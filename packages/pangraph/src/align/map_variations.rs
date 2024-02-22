use crate::align::nextclade::align::params::GapAlignmentSide;
use crate::align::nextclade::align_with_nextclade::{align_with_nextclade, AlignWithNextcladeOutput, NextalignParams};
use crate::align::nextclade::alphabet::nuc::{from_nuc, from_nuc_seq};
use crate::pangraph::edits::{Del, Edits, Ins, Sub};
use eyre::Report;
use itertools::Itertools;

pub fn map_variations(ref_seq: impl AsRef<str>, qry_seq: impl AsRef<str>) -> Result<Edits, Report> {
  let params = NextalignParams {
    min_length: 3,
    ..NextalignParams::default()
  };

  let AlignWithNextcladeOutput {
    substitutions,
    deletions,
    insertions,
    ..
  } = align_with_nextclade(ref_seq, qry_seq, &params)?;

  let subs = substitutions
    .iter()
    .map(|s| Sub::new(s.pos.inner as usize, from_nuc(s.qry_nuc)))
    .collect_vec();

  let dels = deletions
    .iter()
    .map(|d| Del::new(d.range().begin.inner as usize, d.range().end.inner as usize))
    .collect_vec();

  let inss = insertions
    .iter()
    .map(|s| Ins::new(s.pos as usize, from_nuc_seq(&s.ins)))
    .collect_vec();

  Ok(Edits { subs, dels, inss })
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_map_variations_simple_case() {
    // # example alignment
    // #        0            1
    // #        012   3456789012345678
    // # ref = "ACT---TTGCGTATTTACTATA"
    // # qry = "ACTAGATTGAGTATCT---ATA"
    // # sub =           x    x
    // # ins =     xxx
    // # del =                  xxx
    //
    let r = "ACTTTGCGTATTTACTATA";
    let q = "ACTAGATTGAGTATCTATA";

    let actual = map_variations(r, q).unwrap();

    let expected = Edits {
      subs: vec![Sub::new(6, 'A'), Sub::new(11, 'C')],
      dels: vec![Del::new(13, 3)],
      inss: vec![Ins::new(3, "AGA")],
    };

    assert_eq!(expected, actual);

    assert_eq!(r, expected.apply(r).unwrap());
  }
}
