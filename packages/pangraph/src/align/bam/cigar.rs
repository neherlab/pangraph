use eyre::{Report, WrapErr};
use itertools::Itertools;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::cigar::Op;
use noodles::sam::record::Cigar;
use std::str::FromStr;

pub fn parse_cigar_str(cigar_str: impl AsRef<str>) -> Result<Cigar, Report> {
  let cigar_str = cigar_str.as_ref();
  let cigar_str_fixed: String = cigar_str.chars().filter(|c| !"\t ".contains(*c)).collect();
  Cigar::from_str(&cigar_str_fixed).wrap_err_with(|| format!("When parsing CIGAR string '{cigar_str}'"))
}

pub fn cigar_matches_len(cigar: &Cigar) -> usize {
  cigar
    .iter()
    .filter(|op| matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch))
    .map(|op| op.len())
    .sum()
}

pub fn cigar_total_len(cigar: &Cigar) -> usize {
  cigar.iter().map(|op| op.len()).sum()
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use noodles::sam::record::cigar::op::Kind;
  use noodles::sam::record::cigar::Op;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_parse_cigar_string() -> Result<(), Report> {
    let cigar_str = "3H2S10M1I5M1D20M2P10=5X3I4M";

    let actual = parse_cigar_str(cigar_str)?;

    assert_eq!(actual.len(), 12);

    let expected = Cigar::try_from(vec![
      Op::new(Kind::HardClip, 3),
      Op::new(Kind::SoftClip, 2),
      Op::new(Kind::Match, 10),
      Op::new(Kind::Insertion, 1),
      Op::new(Kind::Match, 5),
      Op::new(Kind::Deletion, 1),
      Op::new(Kind::Match, 20),
      Op::new(Kind::Pad, 2),
      Op::new(Kind::SequenceMatch, 10),
      Op::new(Kind::SequenceMismatch, 5),
      Op::new(Kind::Insertion, 3),
      Op::new(Kind::Match, 4),
    ])?;

    assert_eq!(actual, expected);

    Ok(())
  }

  #[rstest]
  fn test_parse_cigar_string_with_tabs() -> Result<(), Report> {
    let cigar_str = "3H\t2S\t10M\t1I\t5M\t1D\t20M\t2P\t10=\t5X\t3I\t4M";

    let actual = parse_cigar_str(cigar_str)?;

    assert_eq!(actual.len(), 12);

    let expected = Cigar::try_from(vec![
      Op::new(Kind::HardClip, 3),
      Op::new(Kind::SoftClip, 2),
      Op::new(Kind::Match, 10),
      Op::new(Kind::Insertion, 1),
      Op::new(Kind::Match, 5),
      Op::new(Kind::Deletion, 1),
      Op::new(Kind::Match, 20),
      Op::new(Kind::Pad, 2),
      Op::new(Kind::SequenceMatch, 10),
      Op::new(Kind::SequenceMismatch, 5),
      Op::new(Kind::Insertion, 3),
      Op::new(Kind::Match, 4),
    ])?;

    assert_eq!(actual, expected);

    Ok(())
  }
}
