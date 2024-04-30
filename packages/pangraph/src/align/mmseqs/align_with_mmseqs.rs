use crate::align::alignment::Alignment;
use crate::align::alignment_args::AlignmentArgs;
use crate::align::mmseqs::paf::PafTsvRecord;
use crate::io::fasta::{write_one_fasta, FastaRecord, FastaWriter};
use crate::io::file::open_file_or_stdin;
use crate::io::fs::path_to_str;
use crate::o;
use crate::utils::subprocess::{create_arg_optional, subprocess, subprocess_with_args};
use cmd_lib::run_cmd;
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::io::Read;
use std::str::FromStr;
use tempfile::Builder as TempDirBuilder;

pub fn align_with_mmseqs(
  refs: &[impl AsRef<str>],
  qrys: &[impl AsRef<str>],
  params: &AlignmentArgs,
) -> Result<Vec<Alignment>, Report> {
  // TODO: This uses a global resource - filesystem.
  // Need to ensure that there are no data races when running concurrently.
  let temp_dir = TempDirBuilder::new().tempdir()?;
  let temp_dir_path = temp_dir.path();
  let qry_path = path_to_str(&temp_dir_path.join("qry.fa"))?;
  let ref_path = path_to_str(&temp_dir_path.join("reff.fa"))?;
  let res_path = path_to_str(&temp_dir_path.join("res.paf"))?;
  let tmp_path = path_to_str(&temp_dir_path.join("tmp"))?;

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

  let output_column_names = PafTsvRecord::fields_names().join(",");
  let k = create_arg_optional("-k", &params.kmer_length);

  run_cmd!(
    mmseqs easy-search
    $qry_path $ref_path $res_path $tmp_path
    --threads 1
    --max-seq-len 10000
    -a
    --search-type 3
    --format-output $output_column_names
    $[k]
  )
  .wrap_err("When trying to align sequences using mmseqs")?;

  let mut paf_str = String::new();
  open_file_or_stdin(&Some(res_path))?.read_to_string(&mut paf_str)?;

  Alignment::from_paf_str(&paf_str)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::pangraph::strand::Strand;
  use crate::utils::interval::Interval;
  use eyre::Report;
  use noodles::sam::record::cigar::op::Kind;
  use noodles::sam::record::cigar::Op;
  use noodles::sam::record::Cigar;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_align_with_mmseqs_general_case() -> Result<(), Report> {
    //                0         1         2         3         4         5         6         7         8         9
    //                0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    let ref_seq = o!("CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTGAGCCTTTGT");
    let qry_seq = o!("CTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATCAAGCTCATGGGACANNNNNCATGGTGGACAGCCTTTGT");

    let params = AlignmentArgs {
      kmer_length: Some(12),
      ..AlignmentArgs::default()
    };

    let actual = align_with_mmseqs(&[ref_seq], &[qry_seq], &params)?;

    let expected = vec![Alignment {
      qry: Hit::new("qry_0", 90, (1, 90)),
      reff: Hit::new("ref_0", 88, (1, 88)),
      matches: 77,
      length: 95,
      quality: 102,
      orientation: Strand::Forward,
      cigar: Cigar::try_from(vec![
        Op::new(Kind::Match, 18),
        Op::new(Kind::Insertion, 4),
        Op::new(Kind::Match, 30),
        Op::new(Kind::Deletion, 5),
        Op::new(Kind::Match, 26),
        Op::new(Kind::Insertion, 3),
        Op::new(Kind::Match, 9),
      ])?,
      divergence: Some(0.18999999999999995),
      align: Some(112.0),
    }];

    assert_eq!(expected, actual);
    Ok(())
  }
}
