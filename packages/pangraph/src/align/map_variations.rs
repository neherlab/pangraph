use crate::align::nextclade::align_with_nextclade::{AlignWithNextcladeOutput, NextalignParams, align_with_nextclade};
use crate::align::nextclade::alphabet::nuc::{from_nuc, from_nuc_seq};
use crate::commands::build::build_args::PangraphBuildArgs;
use crate::pangraph::edits::{Del, Edit, Ins, Sub};
use crate::representation::seq::Seq;
use eyre::Report;
use itertools::Itertools;

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct BandParameters {
  mean_shift: i32,
  band_width: usize,
}

impl BandParameters {
  pub fn new(mean_shift: i32, band_width: usize) -> Self {
    Self { mean_shift, band_width }
  }

  pub fn mean_shift(&self) -> i32 {
    self.mean_shift
  }

  pub fn band_width(&self) -> usize {
    self.band_width
  }

  pub fn add(&mut self, other: &Self) {
    self.mean_shift += other.mean_shift;
    self.band_width += other.band_width;
  }

  pub fn from_edits(edit: &Edit, ref_len: usize) -> Self {
    let mean_shift = edit.aln_mean_shift(ref_len).unwrap();
    let band_width = edit.aln_bandwidth(ref_len, mean_shift).unwrap();
    Self { mean_shift, band_width }
  }
}

pub fn map_variations(
  ref_seq: &Seq,
  qry_seq: &Seq,
  mut band_params: BandParameters,
  args: &PangraphBuildArgs,
) -> Result<Edit, Report> {
  let params = NextalignParams {
    min_length: 1,
    max_alignment_attempts: args.max_alignment_attempts,
    ..NextalignParams::default()
  };

  band_params.band_width += args.extra_band_width;

  let AlignWithNextcladeOutput {
    substitutions,
    deletions,
    insertions,
    ..
  } = align_with_nextclade(ref_seq.as_str(), qry_seq.as_str(), band_params, &params)?;

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
    .map(|s| Ins::new((s.pos + 1) as usize, Seq::from_str(&from_nuc_seq(&s.ins))))
    .collect_vec();

  Ok(Edit { subs, dels, inss })
}

#[cfg(test)]
mod tests {
  use super::*;
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

    let mean_shift = -2;
    let bandwidth = 3;
    let actual = map_variations(
      &Seq::from(r),
      &Seq::from(q),
      BandParameters {
        mean_shift,
        band_width: bandwidth,
      },
      &PangraphBuildArgs::default(),
    )
    .unwrap();

    let expected = Edit {
      subs: vec![Sub::new(6, 'A')],
      dels: vec![Del::new(29, 4)],
      inss: vec![Ins::new(3, "AGA")],
    };

    // test that our example is correct
    let exp_mean_shift = expected.aln_mean_shift(r.len()).unwrap();
    let exp_shift = expected.aln_bandwidth(r.len(), exp_mean_shift).unwrap();
    assert_eq!(exp_mean_shift, mean_shift);
    assert_eq!(exp_shift, bandwidth);

    // test that our example is correct
    assert_eq!(q, &expected.apply(r).unwrap());

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
    let mean_shift = 2;
    let bandwidth = 2;
    let actual = map_variations(
      &Seq::from(r),
      &Seq::from(q),
      BandParameters {
        mean_shift,
        band_width: bandwidth,
      },
      &PangraphBuildArgs::default(),
    )
    .unwrap();

    let expected = Edit {
      subs: vec![Sub::new(10, 'A')],
      dels: vec![Del::new(0, 3), Del::new(36, 4)],
      inss: vec![Ins::new(21, "GGT")],
    };

    // test that our example is correct
    let exp_mean_shift = expected.aln_mean_shift(r.len()).unwrap();
    let exp_shift = expected.aln_bandwidth(r.len(), exp_mean_shift).unwrap();
    assert_eq!(exp_mean_shift, mean_shift);
    assert_eq!(exp_shift, bandwidth);

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
    let mean_shift = -4;
    let bandwidth = 2;
    let actual = map_variations(
      &Seq::from(r),
      &Seq::from(q),
      BandParameters {
        mean_shift,
        band_width: bandwidth,
      },
      &PangraphBuildArgs::default(),
    )
    .unwrap();

    let expected = Edit {
      subs: vec![Sub::new(10, 'A'), Sub::new(31, 'C')],
      dels: vec![Del::new(15, 2)],
      inss: vec![Ins::new(0, "CCTG"), Ins::new(21, "GGT"), Ins::new(40, "GCCGCCG")],
    };

    // test that our example is correct
    let exp_mean_shift = expected.aln_mean_shift(r.len()).unwrap();
    let exp_shift = expected.aln_bandwidth(r.len(), exp_mean_shift).unwrap();
    assert_eq!(exp_mean_shift, mean_shift);
    assert_eq!(exp_shift, bandwidth);

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
    let mean_shift = -2;
    let bandwidth = 11;
    let actual = map_variations(
      &Seq::from(r),
      &Seq::from(q),
      BandParameters {
        mean_shift,
        band_width: bandwidth,
      },
      &PangraphBuildArgs::default(),
    )
    .unwrap();

    let expected = Edit {
      subs: vec![],
      dels: vec![Del::new(21, 9)],
      inss: vec![Ins::new(21, "GGGGGGGGGGGGG")],
    };

    // test that our example is correct
    let exp_mean_shift = expected.aln_mean_shift(r.len()).unwrap();
    let exp_shift = expected.aln_bandwidth(r.len(), exp_mean_shift).unwrap();
    assert_eq!(exp_mean_shift, mean_shift);
    assert_eq!(exp_shift, bandwidth);

    // test that our example is correct
    assert_eq!(q, expected.apply(r).unwrap());

    // test that the aligner reconstructs the variations correctly
    assert_eq!(expected, actual);

    // test that the reconstructed variations are correct
    assert_eq!(q, actual.apply(r).unwrap());
  }
}
