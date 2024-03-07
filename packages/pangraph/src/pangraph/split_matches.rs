use crate::align::alignment::{Alignment, Hit};
use crate::align::alignment_args::AlignmentArgs;
use crate::align::bam::cigar::{cigar_matches_len, cigar_total_len};
use crate::pangraph::strand::Strand;
use crate::{make_internal_error, o};
use eyre::Report;
use itertools::Itertools;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::Cigar;
use serde::{Deserialize, Serialize};
use std::cmp::max;

/// Split the alignments whenever an alignment contains an in/del longer than the threshold length.
pub fn split_matches(aln: &Alignment, args: &AlignmentArgs) -> Result<Vec<Alignment>, Report> {
  let groups = keep_groups(&aln.cigar, args)?;
  let alns = groups
    .into_iter()
    .map(|group| generate_subalignment(aln, &group))
    .collect::<Result<Vec<_>, Report>>()?
    .into_iter()
    .collect_vec();
  Ok(alns)
}

/// Given a cigar string, returns a list of tuples (start_index, end_index) of cigar elements to keep.
/// A kept group must:
///   - always start and end with a match.
///   - have at least `threshold` matches.
///   - have no indels longer than `threshold`.
#[allow(non_snake_case)]
pub fn keep_groups(cigar: &Cigar, args: &AlignmentArgs) -> Result<Vec<(usize, usize)>, Report> {
  let mut groups = vec![];
  let mut g_start = None;
  let mut last_match = None;
  let mut M_sum = 0;
  let mut I_sum = 0;
  let mut D_sum = 0;
  for (i, op) in cigar.iter().enumerate() {
    // Discard leading indels (skip to the first match?)
    if g_start.is_none() {
      if !matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch) {
        continue;
      }
      g_start = Some(i);
    }

    match op.kind() {
      Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
        M_sum += op.len();
        // Reset n. indels from last match
        I_sum = 0;
        D_sum = 0;
        last_match = Some(i);
      }
      Kind::Insertion => {
        I_sum += op.len();
      }
      Kind::Deletion => {
        D_sum += op.len();
      }
      Kind::SoftClip | Kind::HardClip | Kind::Pad | Kind::Skip => {
        // TODO: Should these cases also be handled?
        return make_internal_error!("Unexpected CIGAR operation: '{:#}'", op.kind());
      }
    }

    // Of indel is longer than threshold, split groups
    if max(I_sum, D_sum) >= args.indel_len_threshold {
      if let (Some(g_start), Some(last_match)) = (g_start, last_match) {
        if M_sum >= args.indel_len_threshold {
          // Add only if matches longer than "threshold"
          groups.push((g_start, last_match));
        }
      }
      g_start = None;
      last_match = None;
      M_sum = 0;
      I_sum = 0;
      D_sum = 0;
    }
  }

  // Add last one
  if let (Some(g_start), Some(last_match)) = (g_start, last_match) {
    if M_sum >= args.indel_len_threshold {
      groups.push((g_start, last_match));
    }
  }

  Ok(groups)
}

fn cigar_position_start(cigar: &Cigar, cigar_idx: usize, accepted_ops: &[Kind]) -> usize {
  let mut pos = 0;
  for (i, op) in cigar.iter().enumerate() {
    if i == cigar_idx {
      return pos;
    }
    if accepted_ops.contains(&op.kind()) {
      pos += op.len();
    }
  }
  unreachable!();
}

fn cigar_position_end(cigar: &Cigar, cigar_idx: usize, accepted_ops: &[Kind]) -> usize {
  let mut pos = 0;
  for (i, op) in cigar.iter().enumerate() {
    if accepted_ops.contains(&op.kind()) {
      pos += op.len();
    }
    if i == cigar_idx {
      return pos;
    }
  }
  unreachable!();
}

struct GroupedPositions {
  qry_beg: usize,
  qry_end: usize,
  ref_beg: usize,
  ref_end: usize,
}

fn group_positions(cigar: &Cigar, start_idx: usize, end_idx: usize) -> GroupedPositions {
  let accepted_ops = &[
    Kind::Match,
    Kind::Insertion,
    Kind::SequenceMatch,
    Kind::SequenceMismatch,
  ];
  let qry_beg = cigar_position_start(cigar, start_idx, accepted_ops);
  let qry_end = cigar_position_end(cigar, end_idx, accepted_ops);

  let accepted_ops = &[Kind::Match, Kind::Deletion, Kind::SequenceMatch, Kind::SequenceMismatch];
  let ref_beg = cigar_position_start(cigar, start_idx, accepted_ops);
  let ref_end = cigar_position_end(cigar, end_idx, accepted_ops);

  GroupedPositions {
    qry_beg,
    qry_end,
    ref_beg,
    ref_end,
  }
}

/// Given an alignment and a pair of indices, returns a new alignment that spans only the interval between the two
/// selected cigar indices (python 0-based indexing, right-extreme excluded)
fn generate_subalignment(aln: &Alignment, group: &(usize, usize)) -> Result<Alignment, Report> {
  let GroupedPositions {
    qry_beg: qs,
    qry_end: qe,
    ref_beg: rs,
    ref_end: re,
  } = group_positions(&aln.cigar, group.0, group.1);

  let rs = aln.reff.start + rs;
  let re = aln.reff.start + re;

  let (qs, qe) = match aln.orientation {
    Strand::Forward => (aln.qry.start + qs, aln.qry.start + qe),
    Strand::Reverse => (aln.qry.stop - qe, aln.qry.stop - qs),
  };

  let qry = Hit {
    name: aln.qry.name.clone(),
    length: aln.qry.length,
    start: qs,
    stop: qe,
  };

  let reff = Hit {
    name: aln.reff.name.clone(),
    length: aln.reff.length,
    start: rs,
    stop: re,
  };

  let ops = aln.cigar[group.0..=group.1].iter().copied().collect_vec();
  let cigar = Cigar::try_from(ops)?;

  Ok(Alignment {
    qry,
    reff,
    matches: cigar_matches_len(&cigar),
    length: cigar_total_len(&cigar),
    quality: aln.quality,
    orientation: aln.orientation.clone(),
    cigar,
    divergence: aln.divergence,
    align: aln.align,
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::align::bam::cigar::parse_cigar_str;
  use crate::o;
  use crate::pangraph::strand::Strand;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_keep_groups_simple_case() {
    //        | no                   | keep                 | no      | keep              | no   | keep
    let cig = "10I 20D 10M 20I 190D   40M 1D 1I 40M 1I 40M   1D 100I   200M 60I 60D 140M   200D   40M 2I 70M";
    //           0   1   2   3    4     5  6  7   8  9  10   11   12     13  14  15   16     17    18 19  20

    let expected = vec![(5, 10), (13, 16), (18, 20)];

    let cig = parse_cigar_str(cig.replace(' ', "")).unwrap();
    let args = AlignmentArgs {
      indel_len_threshold: 100,
      ..AlignmentArgs::default()
    };
    let actual = keep_groups(&cig, &args).unwrap();

    assert_eq!(expected, actual);
  }

  #[rstest]
  fn test_split_matches_simple_case_forward() {
    let aln = Alignment {
      qry: Hit {
        name: o!("qry"),
        length: 500,
        start: 200,
        stop: 255,
      },
      reff: Hit {
        name: o!("ref"),
        length: 500,
        start: 100,
        stop: 140,
      },
      matches: 0,
      length: 0,
      quality: 10,
      cigar: parse_cigar_str("3I 6M 3I 3M 4D 5M 14I 7M 3D 4I 5M 5D 3M 3I".replace(' ', "")).unwrap(),
      orientation: Strand::Forward,
      divergence: Some(0.1),
      align: None,
    };

    let actual = split_matches(
      &aln,
      &AlignmentArgs {
        indel_len_threshold: 10,
        ..AlignmentArgs::default()
      },
    )
    .unwrap();

    let expected = vec![
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 203,
          stop: 220,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 100,
          stop: 118,
        },
        matches: 14,
        length: 21,
        quality: 10,
        cigar: parse_cigar_str("6M 3I 3M 4D 5M".replace(' ', "")).unwrap(),
        orientation: Strand::Forward,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 234,
          stop: 253,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 118,
          stop: 141,
        },
        matches: 15,
        length: 27,
        quality: 10,
        cigar: parse_cigar_str("7M 3D 4I 5M 5D 3M".replace(' ', "")).unwrap(),
        orientation: Strand::Forward,
        divergence: Some(0.1),
        align: None,
      },
    ];

    assert_eq!(expected, actual);
  }

  #[rstest]
  fn test_split_matches_simple_case_reverse() {
    let aln = Alignment {
      qry: Hit {
        name: o!("qry"),
        length: 500,
        start: 200,
        stop: 256,
      },
      reff: Hit {
        name: o!("ref"),
        length: 500,
        start: 100,
        stop: 141,
      },
      matches: 0,
      length: 0,
      quality: 10,
      cigar: parse_cigar_str("3I 6M 3I 3M 4D 5M 14I 7M 3D 4I 5M 5D 3M 3I".replace(' ', "")).unwrap(),
      orientation: Strand::Reverse,
      divergence: Some(0.1),
      align: None,
    };

    let actual = split_matches(
      &aln,
      &AlignmentArgs {
        indel_len_threshold: 10,
        ..AlignmentArgs::default()
      },
    )
    .unwrap();

    let expected = vec![
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 236,
          stop: 253,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 100,
          stop: 118,
        },
        matches: 14,
        length: 21,
        quality: 10,
        cigar: parse_cigar_str("6M 3I 3M 4D 5M".replace(' ', "")).unwrap(),
        orientation: Strand::Reverse,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 203,
          stop: 222,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 118,
          stop: 141,
        },
        matches: 15,
        length: 27,
        quality: 10,
        cigar: parse_cigar_str("7M 3D 4I 5M 5D 3M".replace(' ', "")).unwrap(),
        orientation: Strand::Reverse,
        divergence: Some(0.1),
        align: None,
      },
    ];

    assert_eq!(expected, actual);
  }
}
