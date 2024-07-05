use crate::align::alignment::Alignment;
use crate::align::alignment_args::AlignmentArgs;
use crate::align::bam::cigar::parse_cigar_str;
use crate::align::minimap2::minimap2_paf::MinimapPafTsvRecord;
use crate::io::fasta::{write_one_fasta, FastaWriter};
use crate::io::file::{create_file_or_stdout, open_file_or_stdin};
use crate::io::fs::path_to_str;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::utils::subprocess::create_arg_optional;
use cmd_lib::run_cmd;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::Write;
use tempfile::Builder as TempDirBuilder;

pub fn align_with_minimap2(
  blocks: &BTreeMap<BlockId, PangraphBlock>,
  params: &AlignmentArgs,
) -> Result<Vec<Alignment>, Report> {
  // TODO: This uses a global resource - filesystem.
  // Need to ensure that there are no data races when running concurrently.
  let temp_dir = TempDirBuilder::new().tempdir()?;
  let temp_dir_path = temp_dir.path();
  let input_path = path_to_str(&temp_dir_path.join("input.fa"))?;
  let output_path = path_to_str(&temp_dir_path.join("output.paf"))?;

  {
    let mut writer = FastaWriter::from_path(&input_path)?;
    blocks
      .iter()
      .try_for_each(|(id, block)| writer.write(&id.to_string(), block.consensus()))?;
  }

  let kmer_size = create_arg_optional("-k", &params.kmer_length);
  let preset = create_arg_optional("-x", &Some(format!("asm{}", &params.sensitivity)));

  // TODO: implement proper minimap2 call
  #[rustfmt::skip]
  run_cmd!(
    minimap2
    -c
    $input_path
    $input_path
    $[kmer_size]
    $[preset]
    > $output_path
    2> /dev/null
  ).wrap_err("When trying to align sequences using minimap2")?;

  let mut paf_str = String::new();
  open_file_or_stdin(&Some(output_path))?.read_to_string(&mut paf_str)?;

  Alignment::from_minimap_paf_str(&paf_str, blocks)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::o;
  use crate::pangraph::strand::Strand;
  use crate::utils::interval::Interval;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  // TODO: add test cases for alignment of multiple sequences at once

  #[rstest]
  fn test_align_with_minimap2_one_general_case() -> Result<(), Report> {
    let ref_seq = "GTAGTTTTTGTACCCCCCGACTGTACGTCCTCCGTCGATAGAAGCAATAAAGGTGACGTC
    TGACTACTTTTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
    CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTAACATCCGCATCCGGCT
    AACAGCAAGGAAGCCAGGTTTCGGTTCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
    TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCTACTT
    CATCATTTGGATACTTATACTAGCTTGGGCGCTGACTAGTAGTGGCTGCTTGAGCCCCCT
    CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGAGCGCCC
    TACACACCGTCTGATTTCCACATTCCTGCAGGAGCTTACCCGGACCACCGTCAATACACG
    GGATAATAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
    CATCGTTTCATACATGCAAAAGGTTTGGGTGTACCTTATGCTAGGCAGGACCTTTTTAGG
    TGTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATAT
    CGTCTACTATAGGTAAATAGTCCGTTATACATCTCAATTGGAGCGGTTGATGAGCTGACA
    TGTTGGTTACTGTCGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCTAATCATGA
    AAGATAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTG
    CAATCCAGATGACGAAAGCGCTCCCGCCAACGATGACGCAAATAATCGTGCGGTAGGGAG
    TCCTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAG
    GTGGCGCACATTTGATGGTGACTAAGCTGCCAAACTGATT";
    let qry_seq = "GTAGTTCTTGTACCCCCCGACTGTACGTACTCCGTCGATTGAATCAATAAAGGTGACGTC
    TGACTACTATTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
    CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTTACATCCGCATCCGGCT
    AACAGCAAGGAAGCCAGGTTTCCGTGCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
    TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCCACTT
    CATCATTTGGATACTTATACTAGATTGGGCGCTGACTAGTGGTGGCTGCTTGAGCCCCCT
    CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGCGCGCCC
    TACACACCGTCTGACTTCCATATTGCAGCAGGAGCTTACCCGGACCACCGTCAATACACG
    GGATAAGAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
    CATCGCTCATACCTGCAAAAGGTTTGAGTGTACCTTATGCTAGGCAGGGCCTTTTTAGGT
    GTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATATC
    GTCTACCATAGGTAAATAGTCCGTTATACATCTTAATTGGAGCGGTTGATGAGCTGACAT
    GTTGGTTACTGTTGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCAATCATGAAA
    GACAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTGCA
    ATCCAGATGACGAAAGCGCTCCCGCCCACGATGACGCAAATAATCGTGCGGTAGGGAGTC
    CTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAGGT
    GGCGCAGATTTGATGGTGACTAAGCTGCCAAACTGAGT";

    // remove newline characters and spaces
    let refs = vec![ref_seq.replace(['\n', ' '], "")];
    let qrys = vec![qry_seq.replace(['\n', ' '], "")];

    let params = AlignmentArgs {
      kmer_length: Some(10),
      sensitivity: 20,
      ..AlignmentArgs::default()
    };

    let actual = align_with_minimap2(&refs, &qrys, &params)?;

    let expected = vec![Alignment {
      qry: Hit::new("qry_0", 998, (0, 996)),
      reff: Hit::new("ref_0", 1000, (0, 998)),
      matches: 969,
      length: 998,
      quality: 60,
      orientation: Strand::Forward,
      cigar: parse_cigar_str("545M1D225M1D226M").unwrap(),
      divergence: Some(0.0291),
      align: Some(845.0),
    }];

    assert_eq!(expected, actual);
    Ok(())
  }
}
