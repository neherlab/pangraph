use crate::align::alignment::Alignment;
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
  #[clap(long, short = 'd')]
  dummy: Option<usize>,
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
  let tmp_path = path_to_str(&temp_dir_path.join("tmp"))?;

  write_one_fasta(&qry_path, "qry", &qry)?;
  write_one_fasta(&ref_path, "ref", &reff)?;

  let output_column_names = MinimapPafTsvRecord::fields_names().join(",");
  let dummy = create_arg_optional("--dummy", &params.dummy);

  // TODO: implement proper minimap2 call
  #[rustfmt::skip]
  run_cmd!(
    minimap2
    --version
    --query $qry_path
    --ref $ref_path
    $[dummy]
  ).wrap_err("When trying to align sequences using minimap2")?;

  // Write fake data to the output file
  // TODO: remove this once minimap2 call is implemented
  {
    let mut f = create_file_or_stdout(&res_path)?;
    write!(
      f,
      "qry	507	1	497	-	ref	500	500	24	440	508	622	67M10D18M20I235M10I22M1I5M1D119M	0.866	693"
    )?;
  }

  let mut paf_str = String::new();
  open_file_or_stdin(&Some(res_path))?.read_to_string(&mut paf_str)?;

  Alignment::from_minimap_paf_str(&paf_str)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::o;
  use crate::pangraph::pangraph::Strand;
  use eyre::Report;
  use noodles::sam::record::cigar::op::Kind;
  use noodles::sam::record::cigar::Op;
  use noodles::sam::record::Cigar;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_align_with_minimap2_general_case() -> Result<(), Report> {
    //                0         1         2         3         4         5         6         7         8         9
    //                0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    let ref_seq = o!("CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTGAGCCTTTGT");
    let qry_seq = o!("CTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATCAAGCTCATGGGACANNNNNCATGGTGGACAGCCTTTGT");

    let params = Minimap2Params::default();

    let actual = align_with_minimap2(ref_seq, qry_seq, &params)?;

    let expected = Alignment {
      qry: Hit {
        name: o!("qry"),
        length: 90,
        start: 1,
        stop: 90,
        seq: None,
      },
      reff: Hit {
        name: o!("ref"),
        length: 88,
        start: 1,
        stop: 88,
        seq: None,
      },
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
    };

    assert_eq!(expected, actual);
    Ok(())
  }
}
