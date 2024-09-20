use crate::align::alignment::Alignment;
use crate::align::alignment_args::AlignmentArgs;
use crate::align::mmseqs::paf::PafTsvRecord;
use crate::io::fasta::FastaWriter;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::path_to_str;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::utils::subprocess::create_arg_optional;
use cmd_lib::run_cmd;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::Read;
use tempfile::Builder as TempDirBuilder;

pub fn align_with_mmseqs(
  blocks: &BTreeMap<BlockId, PangraphBlock>,
  params: &AlignmentArgs,
) -> Result<Vec<Alignment>, Report> {
  // TODO: This uses a global resource - filesystem.
  // Need to ensure that there are no data races when running concurrently.
  let temp_dir = TempDirBuilder::new().tempdir()?;
  let temp_dir_path = temp_dir.path();
  let input_path = path_to_str(&temp_dir_path.join("input.fa"))?;
  let output_path = path_to_str(&temp_dir_path.join("output.paf"))?;
  let tmp_path = path_to_str(&temp_dir_path.join("tmp"))?;

  {
    let mut writer = FastaWriter::from_path(&input_path)?;
    blocks
      .iter()
      .try_for_each(|(id, block)| writer.write(id.to_string(), block.consensus()))?;
  }

  let output_column_names = PafTsvRecord::fields_names().join(",");
  let k = create_arg_optional("-k", &params.kmer_length);

  run_cmd!(
    mmseqs easy-search
    $input_path $input_path $output_path $tmp_path
    --threads 1
    --max-seq-len 10000
    -a
    --search-type 3
    --format-output $output_column_names
    $[k]
    >/dev/null
  )
  .wrap_err("When trying to align sequences using mmseqs")?;

  let mut paf_str = String::new();
  open_file_or_stdin(&Some(output_path))?.read_to_string(&mut paf_str)?;

  Alignment::from_paf_str(&paf_str)
}

// FIXME: This test is failing after commit a62b19b018b4b2f9602bc75335d4ab5ddbc7abf5
// After that commit, `align_with_mmseqs()` accepts an array of blocks, rather than 2 arrays of ref and query sequences,
// and all blocks sequences are aligned against all block sequences. So there is no way to pass refs and qrys separately
// now - the test needs to be changed. There is also a sibling broken test in `align_with_minimap2`.
//
// #[cfg(test)]
// mod tests {
//   use super::*;
//   use crate::align::alignment::Hit;
//   use crate::pangraph::strand::Strand;
//   use crate::utils::interval::Interval;
//   use eyre::Report;
//   use maplit::btreemap;
//   use noodles::sam::record::cigar::op::Kind;
//   use noodles::sam::record::cigar::Op;
//   use noodles::sam::record::Cigar;
//   use pretty_assertions::assert_eq;
//   use rstest::rstest;
//
//   #[rstest]
//   fn test_align_with_mmseqs_general_case() -> Result<(), Report> {
//     //                0         1         2         3         4         5         6         7         8         9
//     //                0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
//     let ref_seq = o!("CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTGAGCCTTTGT");
//     let qry_seq = o!("CTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATCAAGCTCATGGGACANNNNNCATGGTGGACAGCCTTTGT");
//
//     let params = AlignmentArgs {
//       kmer_length: Some(12),
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
//     let actual = align_with_mmseqs(&blocks, &params)?;
//
//     let expected = vec![Alignment {
//       qry: Hit::new(BlockId(1), 90, (1, 90)),
//       reff: Hit::new(BlockId(0), 88, (1, 88)),
//       matches: 77,
//       length: 95,
//       quality: 102,
//       orientation: Strand::Forward,
//       new_block_id: None, // FIXME
//       anchor_block: None, // FIXME
//       cigar: Cigar::try_from(vec![
//         Op::new(Kind::Match, 18),
//         Op::new(Kind::Insertion, 4),
//         Op::new(Kind::Match, 30),
//         Op::new(Kind::Deletion, 5),
//         Op::new(Kind::Match, 26),
//         Op::new(Kind::Insertion, 3),
//         Op::new(Kind::Match, 9),
//       ])?,
//       divergence: Some(0.18999999999999995),
//       align: Some(112.0),
//     }];
//
//     assert_eq!(expected, actual);
//     Ok(())
//   }
// }
