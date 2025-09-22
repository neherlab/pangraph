use crate::align::map_variations::BandParameters;
use crate::align::nextclade::align::align::align_nuc_simplestripe;
use crate::align::nextclade::align::gap_open::get_gap_open_close_scores_flat;
use crate::align::nextclade::align::insertions_strip::{Insertion, insertions_strip};
pub use crate::align::nextclade::align::params::NextalignParams;
use crate::align::nextclade::alphabet::nuc::{Nuc, from_nuc_seq, to_nuc_seq};
use crate::align::nextclade::analyze::nuc_changes::{FindNucChangesOutput, find_nuc_changes};
use crate::align::nextclade::analyze::nuc_del::NucDelRange;
use crate::align::nextclade::analyze::nuc_sub::NucSub;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct AlignWithNextcladeOutput {
  pub qry_aln: String,
  pub ref_aln: String,
  pub substitutions: Vec<NucSub>,
  pub deletions: Vec<NucDelRange>,
  pub insertions: Vec<Insertion<Nuc>>,
  pub is_reverse_complement: bool,
  pub hit_boundary: bool,
}

pub fn align_with_nextclade(
  reff: impl AsRef<str>,
  qry: impl AsRef<str>,
  band_params: BandParameters,
  params: &NextalignParams,
) -> Result<AlignWithNextcladeOutput, Report> {
  let ref_seq = to_nuc_seq(reff.as_ref()).wrap_err("When converting reference sequence")?;
  let qry_seq = to_nuc_seq(qry.as_ref()).wrap_err("When converting query sequence")?;

  let gap_open_close = get_gap_open_close_scores_flat(&ref_seq, params);

  let alignment = align_nuc_simplestripe(&qry_seq, &ref_seq, &gap_open_close, band_params, params)
    .wrap_err("When aligning sequences with nextclade align")?;

  let stripped = insertions_strip(&alignment.qry_seq, &alignment.ref_seq);

  let FindNucChangesOutput {
    substitutions,
    mut deletions,
    alignment_range,
  } = find_nuc_changes(&stripped.qry_seq, &ref_seq);

  // NB: in nextclade aligner initial/final gaps are not saved as deletions,
  // but they are recorded as limits in the alignment range.
  // We need to add them manually.
  if let Some(alignment_range) = alignment_range {
    if alignment_range.begin.inner > 0 {
      deletions.push(NucDelRange::from_usize(0, alignment_range.begin.inner as usize));
    }
    if (alignment_range.end.inner as usize) < ref_seq.len() {
      deletions.push(NucDelRange::from_usize(
        alignment_range.end.inner as usize,
        ref_seq.len(),
      ));
    }
  } else {
    // If alignment_range is None, it means that there is not interval alignable.
    // We need to add a deletion covering the whole reference sequence.
    deletions.push(NucDelRange::from_usize(0, ref_seq.len()));
  }

  Ok(AlignWithNextcladeOutput {
    qry_aln: from_nuc_seq(&stripped.qry_seq), // returns query with stripped insertions (regions with gap on reference)
    // so that len(qry_aln) == len(ref_seq)
    ref_aln: from_nuc_seq(&alignment.ref_seq), // returns reference sequence with all gaps.
    substitutions,
    deletions,
    insertions: stripped.insertions,
    is_reverse_complement: alignment.is_reverse_complement,
    hit_boundary: alignment.hit_boundary,
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::nextclade::coord::position::NucRefGlobalPosition;
  use crate::commands::build::build_args::PangraphBuildArgs;
  use crate::o;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::sync::LazyLock;

  static EXTRA_BANDWIDTH: LazyLock<usize> = LazyLock::new(|| PangraphBuildArgs::default().extra_band_width);

  #[rstest]
  fn test_align_with_nextclade_general_case() -> Result<(), Report> {
    //                0         1         2         3         4         5         6         7         8         9
    //                0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    let ref_seq = o!("CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTGAGCCTTTGT");
    let qry_seq = o!("CTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATCAAGCTCATGGGACANNTCGCATGGTGGACAGCCTTTGT");
    let ref_aln = o!("CTTGGAGGTTCCGTGGCTA----GATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTG---AGCCTTTGT");
    let qry_aln = o!("CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATC-----AAGCTCATGGGACANNTCGCATGGTGAGCCTTTGT");
    let params = NextalignParams {
      min_length: 3,
      ..NextalignParams::default()
    };

    let band_params = BandParameters::new(0, 4 + *EXTRA_BANDWIDTH);
    let actual = align_with_nextclade(ref_seq, qry_seq, band_params, &params)?;

    let expected = AlignWithNextcladeOutput {
      qry_aln,
      ref_aln,
      substitutions: vec![
        NucSub {
          pos: NucRefGlobalPosition::new(62),
          ref_nuc: Nuc::C,
          qry_nuc: Nuc::G,
        },
        NucSub {
          pos: NucRefGlobalPosition::new(67),
          ref_nuc: Nuc::C,
          qry_nuc: Nuc::N,
        },
        NucSub {
          pos: NucRefGlobalPosition::new(68),
          ref_nuc: Nuc::T,
          qry_nuc: Nuc::N,
        },
      ],
      deletions: vec![NucDelRange::from_usize(48, 53)],
      insertions: vec![
        Insertion {
          pos: 18,
          ins: to_nuc_seq("TAAA")?,
        },
        Insertion {
          pos: 78,
          ins: to_nuc_seq("GAC")?,
        },
      ],
      is_reverse_complement: false,
      hit_boundary: false,
    };

    assert_eq!(expected, actual);

    Ok(())
  }

  #[rstest]
  fn test_align_with_nextclade_ns_in_ref() -> Result<(), Report> {
    //       0         1         2         3         4         5
    //       012345678901234567890123456789012345678901234567890
    // ref = TGGTGCTGCAGCTTATTATGTGGNNNNNTTTTCTATTAAAATATAATGAAA
    // qry = TGGTGCTGCAGCTTATTATGTGGAGGACTTTTCTATTAAAATATAATGAAA
    // sub =                        xxxxx

    let ref_seq = o!("TGGTGCTGCAGCTTATTATGTGGNNNNNTTTTCTATTAAAATATAATGAAA");
    let qry_seq = o!("TGGTGCTGCAGCTTATTATGTGGAGGACTTTTCTATTAAAATATAATGAAA");
    let qry_aln = o!("TGGTGCTGCAGCTTATTATGTGGAGGACTTTTCTATTAAAATATAATGAAA");
    let ref_aln = o!("TGGTGCTGCAGCTTATTATGTGGNNNNNTTTTCTATTAAAATATAATGAAA");

    let params = NextalignParams {
      min_length: 3,
      ..NextalignParams::default()
    };
    let band_params = BandParameters::new(0, *EXTRA_BANDWIDTH);
    let actual = align_with_nextclade(ref_seq, qry_seq, band_params, &params)?;

    let expected = AlignWithNextcladeOutput {
      qry_aln,
      ref_aln,
      substitutions: vec![
        NucSub {
          pos: NucRefGlobalPosition::new(23),
          ref_nuc: Nuc::N,
          qry_nuc: Nuc::A,
        },
        NucSub {
          pos: NucRefGlobalPosition::new(24),
          ref_nuc: Nuc::N,
          qry_nuc: Nuc::G,
        },
        NucSub {
          pos: NucRefGlobalPosition::new(25),
          ref_nuc: Nuc::N,
          qry_nuc: Nuc::G,
        },
        NucSub {
          pos: NucRefGlobalPosition::new(26),
          ref_nuc: Nuc::N,
          qry_nuc: Nuc::A,
        },
        NucSub {
          pos: NucRefGlobalPosition::new(27),
          ref_nuc: Nuc::N,
          qry_nuc: Nuc::C,
        },
      ],
      deletions: vec![],
      insertions: vec![],
      is_reverse_complement: false,
      hit_boundary: false,
    };

    assert_eq!(expected, actual);

    Ok(())
  }

  #[rstest]
  fn test_align_with_nextclade_edge_case() -> Result<(), Report> {
    //       0         1         2         3         4         5         6         7         8         9         0
    //       01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    // ref = TGGTGCTGCNNNNNATTATGTGGGTTATCTTCAACCTTTTTTTAAAATATAATGAAAATGGAACCATTACAGATGCTNNNNNNNNTGCACTTGACCCTCTC
    // qry = TGGTGCTGCAGCTTATTATGTGGGTTATCTTCAACCTTTTTTTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTC
    // sub =          xxxxx                                                               xxxxxxxx

    let ref_seq =
      o!("TGGTGCTGCNNNNNATTATGTGGGTTATCTTCAACCTTTTTTTAAAATATAATGAAAATGGAACCATTACAGATGCTNNNNNNNNTGCACTTGACCCTCTC");
    let qry_seq =
      o!("TGGTGCTGCAGCTTATTATGTGGGTTATCTTCAACCTTTTTTTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTC");

    let qry_aln =
      o!("TGGTGCTGCAGCTTATTATGTGGGTTATCTTCAACCTTTTTTTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTC");
    let ref_aln =
      o!("TGGTGCTGCNNNNNATTATGTGGGTTATCTTCAACCTTTTTTTAAAATATAATGAAAATGGAACCATTACAGATGCTNNNNNNNNTGCACTTGACCCTCTC");

    let params = NextalignParams {
      min_length: 3,
      ..NextalignParams::default()
    };

    let band_params = BandParameters::new(0, *EXTRA_BANDWIDTH);
    let actual = align_with_nextclade(ref_seq, qry_seq, band_params, &params)?;

    let subs = [
      (9, Nuc::A),
      (10, Nuc::G),
      (11, Nuc::C),
      (12, Nuc::T),
      (13, Nuc::T),
      (77, Nuc::G),
      (78, Nuc::T),
      (79, Nuc::A),
      (80, Nuc::G),
      (81, Nuc::A),
      (82, Nuc::C),
      (83, Nuc::T),
      (84, Nuc::G),
    ];

    let subs = subs
      .iter()
      .map(|(pos, nuc)| NucSub {
        pos: NucRefGlobalPosition::new(*pos),
        ref_nuc: Nuc::N,
        qry_nuc: *nuc,
      })
      .collect();

    let expected = AlignWithNextcladeOutput {
      qry_aln,
      ref_aln,
      substitutions: subs,
      deletions: vec![],
      insertions: vec![],
      is_reverse_complement: false,
      hit_boundary: false,
    };

    assert_eq!(expected, actual);

    Ok(())
  }

  #[rstest]
  fn test_align_with_nextclade_unaligned() -> Result<(), Report> {
    // Test case:
    // ref (L = 37) vs qry (L = 18) with bandwidth 0 and mean shift 40
    //       0         1         2         3
    // idx   0123456789012345678901234567890123456
    // ref = AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    // qry = --------------------------------------...----GGGGGGGGGGGGGGGGGG

    let ref_seq = o!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    let qry_seq = o!("GGGGGGGGGGGGGGGGGG");

    let params = NextalignParams {
      min_length: 3,
      ..NextalignParams::default()
    };

    let band_params = BandParameters::new(70, *EXTRA_BANDWIDTH);

    let actual = align_with_nextclade(ref_seq, qry_seq, band_params, &params)?;

    // TODO: What the correct result should be:
    let expected = AlignWithNextcladeOutput {
      qry_aln: o!("-------------------------------------"), // query does not have insertions
      ref_aln: o!("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA------------------"),
      substitutions: vec![],
      deletions: vec![NucDelRange::from_usize(0, 37)], // Should delete entire reference
      insertions: vec![Insertion {
        pos: 36,
        ins: to_nuc_seq("GGGGGGGGGGGGGGGGGG")?,
      }],
      is_reverse_complement: false,
      hit_boundary: false,
    };

    assert_eq!(expected, actual);

    Ok(())
  }
}
