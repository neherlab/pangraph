use crate::align::nextclade::align::params::GapAlignmentSide;
use crate::align::nextclade::align_with_nextclade::{align_with_nextclade, AlignWithNextcladeOutput, NextalignParams};
use crate::align::nextclade::alphabet::nuc::{from_nuc, from_nuc_seq};
use crate::pangraph::edits::{Del, Edit, Ins, Sub};
use eyre::Report;
use itertools::Itertools;

pub fn map_variations(ref_seq: impl AsRef<str>, qry_seq: impl AsRef<str>) -> Result<Edit, Report> {
  let params = NextalignParams {
    min_length: 3,
    min_match_length: 3,
    allowed_mismatches: 16,
    terminal_bandwidth: 100,
    excess_bandwidth: 18,
    max_band_area: 1_000_000_000,
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
    .map(|d| Del::new(d.range().begin.inner as usize, d.range().len()))
    .collect_vec();

  // pangraph convention: location of insertion is the position *after* the insertion -> increment pos by 1
  let inss = insertions
    .iter()
    .map(|s| Ins::new((s.pos + 1) as usize, from_nuc_seq(&s.ins)))
    .collect_vec();

  Ok(Edit { subs, dels, inss })
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_map_variations_simple_case() {
    // example alignment
    //        0            1         2         3
    //        012   345678901234567890123456789012345
    // ref = "ACT---TTGCGTCTGATAGCTTAGCGGATATTGACTGTA"
    // qry = "ACTAGATTGAGTCTGATAGCTTAGCGGATATT----GTA"
    // sub =           x

    let r = "ACTTTGCGTCTGATAGCTTAGCGGATATTTACTGTA";
    let q = "ACTAGATTGAGTCTGATAGCTTAGCGGATATTGTA";

    let actual = map_variations(r, q).unwrap();

    let expected = Edit {
      subs: vec![Sub::new(6, 'A')],
      dels: vec![Del::new(29, 4)],
      inss: vec![Ins::new(3, "AGA")],
    };

    // test that our example is correct
    assert_eq!(q, expected.apply(r).unwrap());

    // test that the aligner reconstructs the variations correctly
    assert_eq!(expected, actual);

    // test that the reconstructed variations are correct
    assert_eq!(q, actual.apply(r).unwrap());
  }

  #[rstest]
  fn test_map_variations_initial_final_deletions() {
    //       0         1         2         3
    //       012345678901234567890   1234567890123456789
    // ref = ACACTGATTTCGTCCCTTAGG---TACTCTACACTGTAGCCTA
    // qry = ---CTGATTTAGTCCCTTAGGGGTTACTCTACACTGTAG----
    // sub =           x

    let r = "ACACTGATTTCGTCCCTTAGGTACTCTACACTGTAGCCTA";
    let q = "CTGATTTAGTCCCTTAGGGGTTACTCTACACTGTAG";

    let actual = map_variations(r, q).unwrap();

    let expected = Edit {
      subs: vec![Sub::new(10, 'A')],
      dels: vec![Del::new(0, 3), Del::new(36, 4)],
      inss: vec![Ins::new(21, "GGT")],
    };

    // test that our example is correct
    assert_eq!(q, expected.apply(r).unwrap());

    // test that the aligner reconstructs the variations correctly
    assert_eq!(expected, actual);

    // test that the reconstructed variations are correct
    assert_eq!(q, actual.apply(r).unwrap());
  }

  #[rstest]
  fn test_map_variations_initial_final_insertions() {
    //           0         1         2            3
    //           012345678901234567890   1234567890123456789
    // ref = ----ACACTGATTTCGTCCCTTAGG---TACTCTACACTGTAGCCTA-------
    // qry = CCTGACACTGATTTAGTCC--TAGGGGTTACTCTACACCGTAGCCTAGCCGCCG
    // sub =               x                       x

    let r = "ACACTGATTTCGTCCCTTAGGTACTCTACACTGTAGCCTA";
    let q = "CCTGACACTGATTTAGTCCTAGGGGTTACTCTACACCGTAGCCTAGCCGCCG";

    let actual = map_variations(r, q).unwrap();

    let expected = Edit {
      subs: vec![Sub::new(10, 'A'), Sub::new(31, 'C')],
      dels: vec![Del::new(15, 2)],
      inss: vec![Ins::new(0, "CCTG"), Ins::new(21, "GGT"), Ins::new(40, "GCCGCCG")],
    };

    // test that our example is correct
    assert_eq!(q, expected.apply(r).unwrap());

    // test that the aligner reconstructs the variations correctly
    assert_eq!(expected, actual);

    // test that the reconstructed variations are correct
    assert_eq!(q, actual.apply(r).unwrap());
  }

  #[rstest]
  fn test_map_variations_overlapping_indels() {
    //       0         1         2                      3         4         5
    //       012345678901234567890             123456789012345678901234567890
    // ref = CGCCCTACTACAAGAGGGAAC-------------TTTTTTTTTAAGTATAGCCACAATAGCTGG
    // qry = CGCCCTACTACAAGAGGGAACGGGGGGGGGGGGG---------AAGTATAGCCACAATAGCTGG

    let r = "CGCCCTACTACAAGAGGGAACTTTTTTTTTAAGTATAGCCACAATAGCTGG";
    let q = "CGCCCTACTACAAGAGGGAACGGGGGGGGGGGGGAAGTATAGCCACAATAGCTGG";

    let actual = map_variations(r, q).unwrap();

    let expected = Edit {
      subs: vec![],
      dels: vec![Del::new(21, 9)],
      inss: vec![Ins::new(21, "GGGGGGGGGGGGG")],
    };

    // test that our example is correct
    assert_eq!(q, expected.apply(r).unwrap());

    // test that the aligner reconstructs the variations correctly
    assert_eq!(expected, actual);

    // test that the reconstructed variations are correct
    assert_eq!(q, actual.apply(r).unwrap());
  }
}
