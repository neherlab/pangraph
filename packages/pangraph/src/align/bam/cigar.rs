use crate::make_error;
use eyre::{Report, WrapErr};
use itertools::Either;
use noodles::sam::record::Cigar;
use noodles::sam::record::cigar::op::{Kind, Op};
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

pub fn invert_cigar(cigar: &Cigar) -> Result<Cigar, Report> {
  let inverted_ops: Vec<Op> = cigar.iter().copied().rev().collect();
  Cigar::try_from(inverted_ops).wrap_err("Failed to create inverted CIGAR")
}

pub fn cigar_switch_ref_qry(cigar: &Cigar) -> Result<Cigar, Report> {
  let switched_ops: Result<Vec<Op>, Report> = cigar
    .iter()
    .map(|op| match op.kind() {
      Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => Ok(*op),
      Kind::Insertion => Ok(Op::new(Kind::Deletion, op.len())),
      Kind::Deletion => Ok(Op::new(Kind::Insertion, op.len())),
      _ => make_error!("CIGAR inversion: unsupported operation kind: {:?}", op.kind()),
    })
    .collect();

  let switched_ops = switched_ops?;

  Cigar::try_from(switched_ops).wrap_err("Failed to create switched CIGAR")
}

pub fn cigar_no_indels(cigar: &Cigar) -> bool {
  cigar
    .iter()
    .all(|op| matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch))
}

pub enum Side {
  Leading,
  Trailing,
}

/// Add a flanking insertion or deletion to the CIGAR.
/// If an in/del is already present before the first match, it will be extended.
/// If there is no in/del before the first match, a new in/del will be added.
/// The side parameter determines whether the in/del is added at the beginning (leading) or end (trailing).
pub fn add_flanking_indel(cigar: &Cigar, kind: Kind, add_len: usize, side: &Side) -> Result<Cigar, Report> {
  if kind != Kind::Insertion && kind != Kind::Deletion {
    return make_error!("Unsupported operation kind for extension: {:?}", kind);
  }

  let cigar_iter = match side {
    Side::Leading => Either::Left(cigar.iter().enumerate()),
    Side::Trailing => Either::Right(cigar.iter().enumerate().rev()),
  };

  let mut replace = None;
  for (i, op) in cigar_iter {
    match op.kind() {
      Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
        break;
      },
      k if k == kind => {
        let new_len = op.len() + add_len;
        let new_op = Op::new(kind, new_len);
        replace = Some((i, new_op));
      },
      _ => (),
    }
  }

  let mut new_ops = cigar.iter().copied().collect::<Vec<Op>>();
  if let Some((i, new_op)) = replace {
    new_ops[i] = new_op;
  } else {
    let ins_pos = match side {
      Side::Leading => 0,
      Side::Trailing => new_ops.len(),
    };
    new_ops.insert(ins_pos, Op::new(kind, add_len));
  }
  Cigar::try_from(new_ops).wrap_err("Failed to create modified CIGAR with flanking indel")
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::utils::error::report_to_string;
  use eyre::Report;
  use noodles::sam::record::cigar::Op;
  use noodles::sam::record::cigar::op::Kind;
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

  #[rstest]
  fn test_invert_cigar() -> Result<(), Report> {
    let cigar_str = "10M1I5M1D20M";
    let cigar = parse_cigar_str(cigar_str)?;

    let inverted = invert_cigar(&cigar)?;

    let expected = Cigar::try_from(vec![
      Op::new(Kind::Match, 20),
      Op::new(Kind::Deletion, 1),
      Op::new(Kind::Match, 5),
      Op::new(Kind::Insertion, 1),
      Op::new(Kind::Match, 10),
    ])?;

    assert_eq!(inverted, expected);

    Ok(())
  }

  #[rstest]
  fn test_switch_ref_qry() -> Result<(), Report> {
    let cigar_str = "10M7I5M1D20M";
    let cigar = parse_cigar_str(cigar_str)?;

    let inverted = cigar_switch_ref_qry(&cigar)?;

    let expected = Cigar::try_from(vec![
      Op::new(Kind::Match, 10),
      Op::new(Kind::Deletion, 7),
      Op::new(Kind::Match, 5),
      Op::new(Kind::Insertion, 1),
      Op::new(Kind::Match, 20),
    ])?;

    assert_eq!(inverted, expected);

    Ok(())
  }

  #[rstest]
  fn test_switch_ref_qry_with_unsupported_op() {
    let cigar_str = "10M2S";
    let cigar = parse_cigar_str(cigar_str).unwrap();

    let result = cigar_switch_ref_qry(&cigar);

    assert_eq!(
      report_to_string(&result.unwrap_err()),
      "CIGAR inversion: unsupported operation kind: SoftClip"
    );
  }

  #[rstest]
  fn test_is_cigar_all_matches() -> Result<(), Report> {
    let cigar_str = "10M20=";
    let cigar = parse_cigar_str(cigar_str)?;

    assert!(cigar_no_indels(&cigar));

    let cigar_str = "10M1I20=";
    let cigar = parse_cigar_str(cigar_str)?;

    assert!(!cigar_no_indels(&cigar));

    Ok(())
  }

  #[rstest]
  fn test_add_leading_indel() -> Result<(), Report> {
    // CIGAR with leading match op
    let cigar = parse_cigar_str("10M5I20M")?;
    let modified = add_flanking_indel(&cigar, Kind::Insertion, 3, &Side::Leading)?;
    // Expect that a new insertion op is added at the beginning since
    // the first op (Match) is not an insertion.
    let expected = Cigar::try_from(vec![
      Op::new(Kind::Insertion, 3),
      Op::new(Kind::Match, 10),
      Op::new(Kind::Insertion, 5),
      Op::new(Kind::Match, 20),
    ])?;
    assert_eq!(modified, expected);

    // Now test when the first op is already the desired indel
    let cigar = parse_cigar_str("5I10M20M")?;
    let modified = add_flanking_indel(&cigar, Kind::Insertion, 3, &Side::Leading)?;
    let expected = Cigar::try_from(vec![
      Op::new(Kind::Insertion, 8),
      Op::new(Kind::Match, 10),
      Op::new(Kind::Match, 20),
    ])?;
    assert_eq!(modified, expected);
    Ok(())
  }

  #[rstest]
  fn test_add_trailing_indel() -> Result<(), Report> {
    // CIGAR with trailing match op
    let cigar = parse_cigar_str("10M5D20M")?;
    let modified = add_flanking_indel(&cigar, Kind::Deletion, 4, &Side::Trailing)?;
    // Expect that a new deletion op is added at the end since
    // the last op (Match) is not a deletion.
    let expected = Cigar::try_from(vec![
      Op::new(Kind::Match, 10),
      Op::new(Kind::Deletion, 5),
      Op::new(Kind::Match, 20),
      Op::new(Kind::Deletion, 4),
    ])?;
    assert_eq!(modified, expected);

    // Now test when the last op is already the desired indel
    let cigar = parse_cigar_str("10M20I")?;
    let modified = add_flanking_indel(&cigar, Kind::Insertion, 4, &Side::Trailing)?;
    let expected = Cigar::try_from(vec![Op::new(Kind::Match, 10), Op::new(Kind::Insertion, 24)])?;
    assert_eq!(modified, expected);
    Ok(())
  }

  #[rstest]
  fn test_add_leading_indel_extend_prefix() -> Result<(), Report> {
    // CIGAR with a prefix having multiple ops before the first match.
    // Here, the prefix is "5D,3I" before the first match, and we extend the insertion.
    let cigar = parse_cigar_str("5D3I10M")?;
    let modified = add_flanking_indel(&cigar, Kind::Insertion, 2, &Side::Leading)?;
    // Expected prefix: "5D" remains, "3I" becomes "5I", then "10M".
    let expected = Cigar::try_from(vec![
      Op::new(Kind::Deletion, 5),
      Op::new(Kind::Insertion, 5),
      Op::new(Kind::Match, 10),
    ])?;
    assert_eq!(modified, expected);
    Ok(())
  }

  #[rstest]
  fn test_add_trailing_indel_extend_suffix() -> Result<(), Report> {
    // CIGAR with trailing region after the last match.
    // For "10M5I3D2I", the head is "10M" and the trailing region is "5I,3D,2I".
    // When adding a trailing deletion, we expect the first deletion op in the suffix to be extended.
    let cigar = parse_cigar_str("10M3D2I")?;
    let modified = add_flanking_indel(&cigar, Kind::Deletion, 4, &Side::Trailing)?;
    // Expected: head "10M", then trailing region becomes "5I, (3D+4=7D),2I".
    let expected = Cigar::try_from(vec![
      Op::new(Kind::Match, 10),
      Op::new(Kind::Deletion, 7),
      Op::new(Kind::Insertion, 2),
    ])?;
    assert_eq!(modified, expected);
    Ok(())
  }

  #[rstest]
  fn test_add_leading_indel_deletion_extend() -> Result<(), Report> {
    // For CIGAR "5D10M", prefix is ["5D"]; when adding deletion type,
    // the prefix should extend to "7D" followed by "10M".
    let cigar = parse_cigar_str("5D10M")?;
    let modified = add_flanking_indel(&cigar, Kind::Deletion, 2, &Side::Leading)?;
    let expected = Cigar::try_from(vec![Op::new(Kind::Deletion, 7), Op::new(Kind::Match, 10)])?;
    assert_eq!(modified, expected);
    Ok(())
  }

  #[rstest]
  fn test_add_trailing_indel_insertion_extend() -> Result<(), Report> {
    // For CIGAR "10M2I", trailing region is ["2I"].
    // When adding an insertion trailing, it should extend "2I" to "5I".
    let cigar = parse_cigar_str("10M2I")?;
    let modified = add_flanking_indel(&cigar, Kind::Insertion, 3, &Side::Trailing)?;
    let expected = Cigar::try_from(vec![Op::new(Kind::Match, 10), Op::new(Kind::Insertion, 5)])?;
    assert_eq!(modified, expected);
    Ok(())
  }
}
