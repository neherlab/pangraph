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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::pangraph::strand::Strand;
  use eyre::Report;
  use maplit::btreemap;
  use minimap2::{Minimap2Index, Minimap2Mapper, Minimap2Options, Minimap2Preset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  // TODO: add test cases for alignment of multiple sequences at once

  #[rstest]
  fn test_align_with_minimap2_lib_one_general_case() -> Result<(), Report> {
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

    // let params = AlignmentArgs {
    //   kmer_length: Some(10),
    //   sensitivity: 20,
    //   ..AlignmentArgs::default()
    // };

    // let blocks = [ref_seq, qry_seq]
    //   .into_iter()
    //   .enumerate()
    //   .map(|(i, seq)| PangraphBlock::new(Some(BlockId(i)), seq.replace(['\n', ' '], ""), btreemap! {}))
    //   .map(|block| (block.id(), block))
    //   .collect();

    let expected = vec![Alignment {
      qry: Hit::new(BlockId(1), 998, (0, 996)),
      reff: Hit::new(BlockId(0), 1000, (0, 998)),
      matches: 969,
      length: 998,
      quality: 60,
      orientation: Strand::Forward,
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("545M1D225M1D226M").unwrap(),
      divergence: Some(0.0291),
      align: Some(845.0),
    }];

    let options = Minimap2Options::with_preset(Minimap2Preset::Asm20)?;
    let idx = Minimap2Index::new(&[ref_seq], &["ref"], options)?;
    let mut mapper = Minimap2Mapper::new(&idx)?;
    let result = mapper.run_map(qry_seq, "qry")?;

    dbg!(&result);

    let actual = parse_cigar_str(&result.regs[0].cigar)?;
    assert_eq!(actual, expected[0].cigar);

    Ok(())
  }
}
