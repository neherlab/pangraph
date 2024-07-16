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
use num_traits::clamp_min;
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

  let min_dp_max = clamp_min(params.indel_len_threshold - 10, 5).to_string();
  let min_dp_max = create_arg_optional("-s", &Some(min_dp_max));

  // TODO: implement proper minimap2 call
  #[rustfmt::skip]
  run_cmd!(
    minimap2
    -c
    -X
    $input_path
    $input_path
    $[kmer_size]
    $[preset]
    $[min_dp_max]
    2> /dev/null
    |
    sed
    "s/\tzd:i:[0-9]*//g"
    |
    sed
    "s/\ts2:i:[0-9]*//g"
    > $output_path
  ).wrap_err("When trying to align sequences using minimap2")?;

  let mut paf_str = String::new();
  open_file_or_stdin(&Some(output_path))?.read_to_string(&mut paf_str)?;

  Alignment::from_minimap_paf_str(&paf_str)
}

// FIXME: This test is failing after commit a62b19b018b4b2f9602bc75335d4ab5ddbc7abf5
// After that commit, `align_with_minimap2()` accepts an array of blocks, rather than 2 arrays of ref and query sequences,
// and all blocks sequences are aligned against all block sequences. So there is no way to pass refs and qrys separately
// now - the test needs to be changed. There is also a sibling broken test in `align_with_mmseqs`.
//
// #[cfg(test)]
// mod tests {
//   use super::*;
//   use crate::align::alignment::Hit;
//   use crate::pangraph::strand::Strand;
//   use eyre::Report;
//   use maplit::btreemap;
//   use pretty_assertions::assert_eq;
//   use rstest::rstest;
//
//   // TODO: add test cases for alignment of multiple sequences at once
//
//   #[rstest]
//   fn test_align_with_minimap2_one_general_case() -> Result<(), Report> {
//     let ref_seq = "GTAGTTTTTGTACCCCCCGACTGTACGTCCTCCGTCGATAGAAGCAATAAAGGTGACGTC
//     TGACTACTTTTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
//     CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTAACATCCGCATCCGGCT
//     AACAGCAAGGAAGCCAGGTTTCGGTTCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
//     TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCTACTT
//     CATCATTTGGATACTTATACTAGCTTGGGCGCTGACTAGTAGTGGCTGCTTGAGCCCCCT
//     CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGAGCGCCC
//     TACACACCGTCTGATTTCCACATTCCTGCAGGAGCTTACCCGGACCACCGTCAATACACG
//     GGATAATAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
//     CATCGTTTCATACATGCAAAAGGTTTGGGTGTACCTTATGCTAGGCAGGACCTTTTTAGG
//     TGTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATAT
//     CGTCTACTATAGGTAAATAGTCCGTTATACATCTCAATTGGAGCGGTTGATGAGCTGACA
//     TGTTGGTTACTGTCGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCTAATCATGA
//     AAGATAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTG
//     CAATCCAGATGACGAAAGCGCTCCCGCCAACGATGACGCAAATAATCGTGCGGTAGGGAG
//     TCCTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAG
//     GTGGCGCACATTTGATGGTGACTAAGCTGCCAAACTGATT";
//     let qry_seq = "GTAGTTCTTGTACCCCCCGACTGTACGTACTCCGTCGATTGAATCAATAAAGGTGACGTC
//     TGACTACTATTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
//     CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTTACATCCGCATCCGGCT
//     AACAGCAAGGAAGCCAGGTTTCCGTGCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
//     TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCCACTT
//     CATCATTTGGATACTTATACTAGATTGGGCGCTGACTAGTGGTGGCTGCTTGAGCCCCCT
//     CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGCGCGCCC
//     TACACACCGTCTGACTTCCATATTGCAGCAGGAGCTTACCCGGACCACCGTCAATACACG
//     GGATAAGAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
//     CATCGCTCATACCTGCAAAAGGTTTGAGTGTACCTTATGCTAGGCAGGGCCTTTTTAGGT
//     GTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATATC
//     GTCTACCATAGGTAAATAGTCCGTTATACATCTTAATTGGAGCGGTTGATGAGCTGACAT
//     GTTGGTTACTGTTGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCAATCATGAAA
//     GACAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTGCA
//     ATCCAGATGACGAAAGCGCTCCCGCCCACGATGACGCAAATAATCGTGCGGTAGGGAGTC
//     CTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAGGT
//     GGCGCAGATTTGATGGTGACTAAGCTGCCAAACTGAGT";
//
//     let params = AlignmentArgs {
//       kmer_length: Some(10),
//       sensitivity: 20,
//       ..AlignmentArgs::default()
//     };
//
//     let blocks = [ref_seq, qry_seq]
//       .into_iter()
//       .enumerate()
//       .map(|(i, seq)| PangraphBlock::new(Some(BlockId(i)), seq.replace(['\n', ' '], ""), btreemap! {}))
//       .map(|block| (block.id(), block))
//       .collect();
//
//     let actual = align_with_minimap2(&blocks, &params)?;
//
//     let expected = vec![Alignment {
//       qry: Hit::new(BlockId(1), 998, (0, 996)),
//       reff: Hit::new(BlockId(0), 1000, (0, 998)),
//       matches: 969,
//       length: 998,
//       quality: 60,
//       orientation: Strand::Forward,
//       new_block_id: None, // FIXME
//       anchor_block: None, // FIXME
//       cigar: parse_cigar_str("545M1D225M1D226M").unwrap(),
//       divergence: Some(0.0291),
//       align: Some(845.0),
//     }];
//
//     assert_eq!(expected, actual);
//     Ok(())
//   }
// }
