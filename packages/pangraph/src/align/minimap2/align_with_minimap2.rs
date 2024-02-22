use crate::align::alignment::Alignment;
use crate::align::bam::cigar::parse_cigar_str;
use crate::align::minimap2::minimap2_paf::MinimapPafTsvRecord;
use crate::io::fasta::write_one_fasta;
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
  kmersize: Option<usize>,
  #[clap(long, short = 'x')]
  preset: Option<String>,
}

pub fn align_with_minimap2(
  reff: impl AsRef<str>,
  qry: impl AsRef<str>,
  params: &Minimap2Params,
) -> Result<Alignment, Report> {
  // TODO: This uses a global resource - filesystem.
  // Need to ensure that there are no data races when running concurrently.
  let temp_dir = TempDirBuilder::new().tempdir()?;
  let temp_dir_path = temp_dir.path();
  let qry_path = path_to_str(&temp_dir_path.join("qry.fa"))?;
  let ref_path = path_to_str(&temp_dir_path.join("reff.fa"))?;
  let res_path = path_to_str(&temp_dir_path.join("res.paf"))?;

  write_one_fasta(&qry_path, "qry", &qry)?;
  write_one_fasta(&ref_path, "ref", &reff)?;

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
  use noodles::sam::record::cigar::op::Kind;
  use noodles::sam::record::cigar::Op;
  use noodles::sam::record::Cigar;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_align_with_minimap2_general_case() -> Result<(), Report> {
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
    let ref_seq = ref_seq.replace("\n", "").replace(" ", "");
    let qry_seq = qry_seq.replace("\n", "").replace(" ", "");

    let params = Minimap2Params {
      kmersize: Some(10),
      preset: Some("asm20".to_string()),
    };

    let actual = align_with_minimap2(ref_seq, qry_seq, &params)?;

    let expected = Alignment {
      qry: Hit {
        name: o!("qry"),
        length: 998,
        start: 0,
        stop: 996,
        seq: None,
      },
      reff: Hit {
        name: o!("ref"),
        length: 1000,
        start: 0,
        stop: 998,
        seq: None,
      },
      matches: 969,
      length: 998,
      quality: 60,
      orientation: Strand::Forward,
      cigar: parse_cigar_str("545M1D225M1D226M").unwrap(),
      divergence: Some(0.0291),
      align: Some(845.0),
    };

    assert_eq!(expected, actual);
    Ok(())
  }
}
