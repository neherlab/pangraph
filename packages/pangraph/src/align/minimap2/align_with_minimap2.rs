use crate::align::alignment::Alignment;
use crate::align::bam::cigar::parse_cigar_str;
use crate::align::minimap2::minimap2_paf::MinimapPafTsvRecord;
use crate::io::fasta::{write_one_fasta, FastaWriter};
use crate::io::file::{create_file_or_stdout, open_file_or_stdin};
use crate::io::fs::path_to_str;
use crate::utils::subprocess::create_arg_optional;
use clap::Args;
use cmd_lib::run_cmd;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use std::io::Write;
use tempfile::Builder as TempDirBuilder;

#[derive(Clone, Debug, Default, Args, Serialize, Deserialize)]
pub struct Minimap2Params {
  #[clap(long, short = 'k')]
  pub kmersize: Option<usize>,
  #[clap(long, short = 'x')]
  pub preset: Option<String>,
}

pub fn align_with_minimap2(
  refs: &[impl AsRef<str>],
  qrys: &[impl AsRef<str>],
  params: &Minimap2Params,
) -> Result<Vec<Alignment>, Report> {
  // TODO: This uses a global resource - filesystem.
  // Need to ensure that there are no data races when running concurrently.
  let temp_dir = TempDirBuilder::new().tempdir()?;
  let temp_dir_path = temp_dir.path();
  let qry_path = path_to_str(&temp_dir_path.join("qry.fa"))?;
  let ref_path = path_to_str(&temp_dir_path.join("reff.fa"))?;
  let res_path = path_to_str(&temp_dir_path.join("res.paf"))?;

  {
    let mut writer = FastaWriter::from_path(&qry_path)?;
    qrys
      .iter()
      .enumerate()
      .try_for_each(|(i, seq)| writer.write(format!("qry_{i}"), seq))?;
  }

  {
    let mut writer = FastaWriter::from_path(&ref_path)?;
    refs
      .iter()
      .enumerate()
      .try_for_each(|(i, seq)| writer.write(format!("ref_{i}"), seq))?;
  }

  let output_column_names = MinimapPafTsvRecord::fields_names().join(",");
  let kmer_size = create_arg_optional("-k", &params.kmersize);
  let preset = create_arg_optional("-x", &params.preset);

  // TODO: implement proper minimap2 call
  #[rustfmt::skip]
  run_cmd!(
    minimap2
    -c
    $ref_path
    $qry_path
    $[kmer_size]
    $[preset]
    > $res_path
    2> /dev/null
  ).wrap_err("When trying to align sequences using minimap2")?;

  let mut paf_str = String::new();
  open_file_or_stdin(&Some(res_path))?.read_to_string(&mut paf_str)?;

  Alignment::from_minimap_paf_str(&paf_str)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::o;
  use crate::pangraph::strand::Strand;
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

    let params = Minimap2Params {
      kmersize: Some(10),
      preset: Some(o!("asm20")),
    };

    let actual = align_with_minimap2(&refs, &qrys, &params)?;

    let expected = vec![Alignment {
      qry: Hit {
        name: o!("qry_0"),
        length: 998,
        start: 0,
        stop: 996,
      },
      reff: Hit {
        name: o!("ref_0"),
        length: 1000,
        start: 0,
        stop: 998,
      },
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
