use crate::align::alignment::{Alignment, Hit};
use crate::align::alignment_args::AlignmentArgs;
use crate::align::bam::cigar::{cigar_matches_len, cigar_total_len};
use crate::make_internal_error;
use crate::pangraph::strand::Strand;
use eyre::Report;
use itertools::Itertools;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::cigar::Op;
use noodles::sam::record::Cigar;
use serde::{Deserialize, Serialize};
use std::cmp::max;

/// Split the alignments whenever an alignment contains an in/del longer than the threshold length.
pub fn split_matches(aln: &Alignment, args: &AlignmentArgs) -> Result<Vec<Alignment>, Report> {
  let groups = keep_groups(&aln.cigar, args)?;

  let mut alns = groups
    .into_iter()
    .map(|group| generate_subalignment(aln, &group))
    .collect::<Result<Vec<Alignment>, Report>>()?;

  alns.iter_mut().try_for_each(|aln| side_patches(aln, args))?;

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

  let rs = aln.reff.interval.start + rs;
  let re = aln.reff.interval.start + re;

  let (qs, qe) = match aln.orientation {
    Strand::Forward => (aln.qry.interval.start + qs, aln.qry.interval.start + qe),
    Strand::Reverse => (aln.qry.interval.end - qe, aln.qry.interval.end - qs),
  };

  let qry = Hit::new(aln.qry.name.clone(), aln.qry.length, (qs, qe));
  let reff = Hit::new(aln.reff.name.clone(), aln.reff.length, (rs, re));
  let ops = aln.cigar[group.0..=group.1].iter().copied().collect_vec();
  let cigar = Cigar::try_from(ops)?;

  Ok(Alignment {
    qry,
    reff,
    matches: cigar_matches_len(&cigar),
    length: cigar_total_len(&cigar),
    quality: aln.quality,
    orientation: aln.orientation,
    cigar,
    divergence: aln.divergence,
    align: aln.align,
  })
}

/// If lateral overhangs are shorter than threshold, add them to the alignment to avoid excessive fragmentation
#[allow(non_snake_case)]
pub fn side_patches(aln: &mut Alignment, args: &AlignmentArgs) -> Result<(), Report> {
  let mut ops = aln.cigar.to_vec();

  // Check reference
  let (rs, re, rL) = (aln.reff.interval.start, aln.reff.interval.end, aln.reff.length);
  if (rs > 0) && (rs < args.indel_len_threshold) {
    // Append left side reference patch
    let delta_l = rs;
    aln.reff.interval.start = 0;
    aln.length += delta_l;
    ops.insert(0, Op::new(Kind::Deletion, delta_l));
  }
  if (re < rL) && (rL - re < args.indel_len_threshold) {
    // Append right reference patch
    let delta_l = rL - re;
    aln.reff.interval.end = rL;
    aln.length += delta_l;
    ops.push(Op::new(Kind::Deletion, delta_l));
  }

  // Check query
  let (qs, qe, qL) = (aln.qry.interval.start, aln.qry.interval.end, aln.qry.length);
  if (qs > 0) && (qs < args.indel_len_threshold) {
    // Append query start
    let delta_l = qs;
    aln.qry.interval.start = 0;
    aln.length += delta_l;
    let extra_ins = Op::new(Kind::Insertion, delta_l);
    if aln.orientation == Strand::Forward {
      ops.push(extra_ins);
    } else {
      ops.insert(0, extra_ins);
    }
  }
  if (qe < qL) && (qL - qe < args.indel_len_threshold) {
    // Append query end
    let delta_l = qL - qe;
    aln.qry.interval.end = qL;
    aln.length += delta_l;
    let extra_ins = Op::new(Kind::Insertion, delta_l);
    if aln.orientation == Strand::Forward {
      ops.push(extra_ins);
    } else {
      ops.insert(0, extra_ins);
    }
  }

  aln.cigar = Cigar::try_from(ops)?;

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::align::bam::cigar::parse_cigar_str;
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
    // CG      3I  6M     3I  3M  4D   5M    14I            7M      3D  4I   5M    5D    3M  3I
    //
    // 100 +       0                      17                18                            40
    // ref     --- MMMMMM --- MMM DDDD MMMMM -------------- MMMMMMM DDD ---- MMMMM DDDDD MMM ---
    // qry     III MMMMMM III MMM ---- MMMMM IIIIIIIIIIIIII MMMMMMM --- IIII MMMMM ----- MMM III
    // 200 +       3                      19                34                            52
    // groups      |-----------------------|                |------------------------------|

    let aln = Alignment {
      qry: Hit::new("qry", 500, (200, 255)),
      reff: Hit::new("ref", 500, (100, 140)),
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
        qry: Hit::new("qry", 500, (203, 220)),
        reff: Hit::new("ref", 500, (100, 118)),
        matches: 14,
        length: 21,
        quality: 10,
        cigar: parse_cigar_str("6M 3I 3M 4D 5M".replace(' ', "")).unwrap(),
        orientation: Strand::Forward,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit::new("qry", 500, (234, 253)),
        reff: Hit::new("ref", 500, (118, 141)),
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
    // CG      3I  6M     3I  3M  4D   5M    14I            7M      3D  4I   5M    5D    3M  3I
    //
    // 100 +       0                      17                18                            40
    // ref     --- MMMMMM --- MMM DDDD MMMMM -------------- MMMMMMM DDD ---- MMMMM DDDDD MMM ---
    // qry     III MMMMMM III MMM ---- MMMMM IIIIIIIIIIIIII MMMMMMM --- IIII MMMMM ----- MMM III
    // 200 +       52                     36                21                             3
    // groups      |-----------------------|                |------------------------------|

    let aln = Alignment {
      qry: Hit::new("qry", 500, (200, 256)),
      reff: Hit::new("ref", 500, (100, 141)),
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
        qry: Hit::new("qry", 500, (236, 253)),
        reff: Hit::new("ref", 500, (100, 118)),
        matches: 14,
        length: 21,
        quality: 10,
        cigar: parse_cigar_str("6M 3I 3M 4D 5M".replace(' ', "")).unwrap(),
        orientation: Strand::Reverse,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit::new("qry", 500, (203, 222)),
        reff: Hit::new("ref", 500, (118, 141)),
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

  #[rstest]
  fn test_split_matches_with_side_patches_forward() {
    // CG      3I  3D  6M     3I  3M  4D   5M    14I            7M      3D  4I   5M    5D    3M  4I   12D
    //
    //             0   3                      20                21                            43                55
    // ref     --- DDD MMMMMM --- MMM DDDD MMMMM -------------- MMMMMMM DDD ---- MMMMM DDDDD MMM ---- DDDDDDDDDDDD
    // qry     III --- MMMMMM III MMM ---- MMMMM IIIIIIIIIIIIII MMMMMMM --- IIII MMMMM ----- MMM IIII ------------
    // 200 +           3                      19                34                            52   56
    // groups          |-----------------------|                |------------------------------|
    // side patch  |---------------------------|                |-----------------------------------|

    let aln = Alignment {
      qry: Hit::new("qry", 257, (200, 257)),
      reff: Hit::new("ref", 56, (0, 56)),
      matches: 29,
      length: 84,
      quality: 10,
      cigar: parse_cigar_str("3I 3D 6M 3I 3M 4D 5M 14I 7M 3D 4I 5M 5D 3M 4I 12D".replace(' ', "")).unwrap(),
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
        qry: Hit::new("qry", 257, (203, 220)),
        reff: Hit::new("ref", 56, (0, 21)),
        matches: 14,
        length: 24,
        quality: 10,
        cigar: parse_cigar_str("3D 6M 3I 3M 4D 5M".replace(' ', "")).unwrap(),
        orientation: Strand::Forward,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit::new("qry", 257, (234, 257)),
        reff: Hit::new("ref", 56, (21, 44)),
        matches: 15,
        length: 31,
        quality: 10,
        cigar: parse_cigar_str("7M 3D 4I 5M 5D 3M 4I".replace(' ', "")).unwrap(),
        orientation: Strand::Forward,
        divergence: Some(0.1),
        align: None,
      },
    ];

    assert_eq!(expected, actual);
  }

  #[rstest]
  fn test_split_matches_with_side_patches_reverse() {
    // CG          3I  3D  6M     3I  3M  4D   5M    14I            7M      3D  4I   5M    5D    3M  4I   5D
    //
    //                 0   3                      20                21                            43         48
    // ref         --- DDD MMMMMM --- MMM DDDD MMMMM -------------- MMMMMMM DDD ---- MMMMM DDDDD MMM ---- DDDDD
    // qry         III --- MMMMMM III MMM ---- MMMMM IIIIIIIIIIIIII MMMMMMM --- IIII MMMMM ----- MMM IIII -----
    // 300 +       56      53                     37                22                             4    0
    // groups              |-----------------------|                |------------------------------|
    // side patch  |-------------------------------|                |-------------------------------      -----|

    let aln = Alignment {
      qry: Hit::new("qry", 257, (200, 257)),
      reff: Hit::new("ref", 49, (0, 49)),
      matches: 29,
      length: 77,
      quality: 10,
      cigar: parse_cigar_str("3I 3D 6M 3I 3M 4D 5M 14I 7M 3D 4I 5M 5D 3M 4I 5D".replace(' ', "")).unwrap(),
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
        qry: Hit::new("qry", 257, (237, 257)),
        reff: Hit::new("ref", 49, (0, 21)),
        matches: 14,
        length: 27,
        quality: 10,
        cigar: parse_cigar_str("3I 3D 6M 3I 3M 4D 5M".replace(' ', "")).unwrap(),
        orientation: Strand::Reverse,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit::new("qry", 257, (204, 223)),
        reff: Hit::new("ref", 49, (21, 49)),
        matches: 15,
        length: 32,
        quality: 10,
        cigar: parse_cigar_str("7M 3D 4I 5M 5D 3M 5D".replace(' ', "")).unwrap(),
        orientation: Strand::Reverse,
        divergence: Some(0.1),
        align: None,
      },
    ];

    assert_eq!(expected, actual);
  }
}
