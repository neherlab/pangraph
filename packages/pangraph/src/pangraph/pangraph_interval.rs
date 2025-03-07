use crate::align::alignment::ExtractedHit;
use crate::make_internal_error;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::strand::Strand;
use crate::utils::id::id;
use crate::utils::interval::Interval;
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use noodles::sam::record::Cigar;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PangraphInterval {
  pub interval: Interval,
  pub aligned: bool,
  pub new_block_id: BlockId,
  pub is_anchor: Option<bool>,
  pub orientation: Option<Strand>,
  pub cigar: Option<Cigar>,
}

impl PangraphInterval {
  pub fn len(&self) -> usize {
    self.interval.len()
  }

  pub fn is_empty(&self) -> bool {
    self.interval.is_empty()
  }

  pub fn contains(&self, pos: usize) -> bool {
    self.interval.contains(pos)
  }

  pub fn has_overlap_with(&self, other: &Interval) -> bool {
    self.interval.has_overlap_with(other)
  }

  pub fn insertion_overlap(&self, ins_pos: usize, block_len: usize) -> bool {
    self.interval.contains(ins_pos) || (ins_pos == block_len && self.interval.end == block_len)
  }
}

pub fn have_no_overlap(intervals: &[PangraphInterval], candidate: &PangraphInterval) -> bool {
  !intervals
    .iter()
    .any(|interval| interval.has_overlap_with(&candidate.interval))
}

#[allow(clippy::missing_asserts_for_indexing)]
fn intervals_sanity_checks(intervals: &[PangraphInterval], block_length: usize) -> Result<(), Report> {
  if intervals.is_empty() {
    return make_internal_error!("Intervals array cannot be empty.");
  }

  if intervals[0].interval.start != 0 {
    return make_internal_error!("First interval does not start at 0: {}", intervals[0].interval);
  }

  if intervals.last().unwrap().interval.end != block_length {
    return make_internal_error!(
      "Last interval does not end at the block length. Interval: {}, block length: {}",
      intervals.last().unwrap().interval,
      block_length
    );
  }

  for n in 1..intervals.len() {
    if intervals[n - 1].interval.end != intervals[n].interval.start {
      return make_internal_error!(
        "Intervals {} and {} are not contiguous: {} and {}",
        n - 1,
        n,
        intervals[n - 1].interval,
        intervals[n].interval
      );
    }

    if !intervals[n - 1].aligned && !intervals[n].aligned {
      return make_internal_error!(
        "Found two consecutive unaligned intervals: {} and {}: {} and {}",
        n - 1,
        n,
        intervals[n - 1].aligned,
        intervals[n].aligned
      );
    }
  }
  Ok(())
}

fn unaligned_interval(interval: Interval, block_id: BlockId) -> PangraphInterval {
  // FIXME: partially initialized object in incorrect/unknown state.
  // We probably don't want to calculate this hash ad-hoc. Refactor code to avoid this.
  let new_block_id = calculate_hash(block_id, &interval);
  PangraphInterval {
    interval,
    aligned: false,
    new_block_id,
    is_anchor: None,
    orientation: None,
    cigar: None,
  }
}

fn aligned_interval(h: &ExtractedHit) -> PangraphInterval {
  // FIXME: partially initialized object in incorrect/unknown state.
  // We probably don't want to calculate this hash ad-hoc. Refactor code to avoid this.
  PangraphInterval {
    interval: h.hit.interval.clone(),
    aligned: true,
    new_block_id: h.new_block_id,
    is_anchor: Some(h.is_anchor),
    orientation: Some(h.orientation),
    cigar: h.cigar.clone(),
  }
}

// FIXME: this should not exist
fn calculate_hash(block_id: BlockId, interval: &Interval) -> BlockId {
  BlockId(id((block_id, interval)))
}

fn create_intervals(hits: &[ExtractedHit], block_length: usize) -> Vec<PangraphInterval> {
  let mut intervals = vec![];
  let mut cursor = 0;
  for h in hits.iter().sorted_by_key(|x| x.hit.interval.start) {
    let hit = &h.hit;
    if hit.interval.start > cursor {
      intervals.push(unaligned_interval(Interval::new(cursor, hit.interval.start), hit.name));
      // cursor = hit.interval.start; // FIXME: value assigned but never read
    }
    intervals.push(aligned_interval(h));
    cursor = hit.interval.end;
  }

  if cursor < block_length {
    intervals.push(unaligned_interval(
      Interval::new(cursor, block_length),
      hits.last().unwrap().hit.name,
    ));
  }

  intervals
}

fn check_interval_conditions(
  interval: &PangraphInterval,
  n: usize,
  left_len: usize,
  right_len: usize,
  thr_len: usize,
  intervals: &[PangraphInterval],
) -> Result<(), Report> {
  if interval.aligned {
    return make_internal_error!("Aligned interval at index {n} is shorter than the threshold length ({thr_len}).")
      .with_section(|| format!("{interval:#?}").header("Interval:"));
  }

  if n > 0 {
    if !intervals[n - 1].aligned {
      return make_internal_error!("No adjacent aligned interval on the left of index {n}.")
        .with_section(|| format!("{:#?}", intervals[n - 1]).header("Interval:"));
    }
    if left_len < thr_len {
      return make_internal_error!(
        "Left interval at index {} is shorter than the threshold length ({thr_len}).",
        n - 1
      )
      .with_section(|| format!("{:#?}", intervals[n - 1]).header("Interval:"));
    }
  }
  if n + 1 < intervals.len() {
    if !intervals[n + 1].aligned {
      return make_internal_error!("No adjacent aligned interval on the right of index {n}.")
        .with_section(|| format!("{:#?}", intervals[n + 1]).header("Interval:"));
    }
    if right_len < thr_len {
      return make_internal_error!(
        "Right interval at index {} is shorter than the threshold length ({thr_len}).",
        n + 1
      )
      .with_section(|| format!("{:#?}", intervals[n + 1]).header("Interval:"));
    }
  }
  Ok(())
}

fn refine_intervals(intervals: &mut Vec<PangraphInterval>, thr_len: usize) -> Result<(), Report> {
  let mut mergers = vec![];

  for (n, interval) in intervals.iter().enumerate() {
    if interval.len() < thr_len {
      let left_len = if n > 0 { intervals[n - 1].len() } else { 0 };
      let right_len = if n + 1 < intervals.len() {
        intervals[n + 1].len()
      } else {
        0
      };

      check_interval_conditions(interval, n, left_len, right_len, thr_len, intervals)?;

      mergers.push(if left_len >= right_len { (n, n - 1) } else { (n, n + 1) });
    }
  }

  for (n_from, n_to) in mergers.iter().rev() {
    if n_from < n_to {
      intervals[*n_to].interval.start = intervals[*n_from].interval.start;
    } else {
      intervals[*n_to].interval.end = intervals[*n_from].interval.end;
    }
    intervals.remove(*n_from);
  }

  Ok(())
}

/// given a list of hits on a block, split the block into a list of intervals.
/// Each interval is either aligned to a new block or unaligned.
/// The intervals are refined by merging unaligned intervals shorter than a threshold length.
pub fn extract_intervals(
  hits: &[ExtractedHit],
  block_length: usize,
  thr_len: usize,
) -> Result<Vec<PangraphInterval>, Report> {
  let mut intervals = create_intervals(hits, block_length);
  refine_intervals(&mut intervals, thr_len)?;
  intervals_sanity_checks(&intervals, block_length)
    .wrap_err("When performing interval sanity checks")
    .with_section(|| format!("{intervals:#?}").header("Intervals:"))?;
  Ok(intervals)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::pangraph::strand::Strand::Forward;
  use pretty_assertions::assert_eq;

  fn example() -> (Vec<ExtractedHit>, usize, BlockId) {
    let block_length = 1000;
    let bid = BlockId(0);

    let create_hit = |new_bid: BlockId, is_anchor: bool, strand: Strand, interval: Interval| -> ExtractedHit {
      ExtractedHit {
        new_block_id: new_bid,
        is_anchor,
        orientation: strand,
        hit: Hit {
          name: bid,
          interval,
          length: 0,
        },
        cigar: None,
      }
    };

    let hits = vec![
      create_hit(BlockId(1), true, Forward, Interval::new(10, 100)),
      create_hit(BlockId(2), false, Forward, Interval::new(200, 300)),
      create_hit(BlockId(3), true, Forward, Interval::new(310, 500)),
      create_hit(BlockId(4), false, Forward, Interval::new(600, 900)),
    ];

    (hits, block_length, bid)
  }

  #[test]
  fn test_create_intervals() {
    let (hits, block_length, bid) = example();
    let intervals = create_intervals(&hits, block_length);

    assert_eq!(
      intervals,
      vec![
        PangraphInterval {
          interval: Interval::new(0, 10),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(0, 10)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(10, 100),
          aligned: true,
          new_block_id: BlockId(1),
          is_anchor: Some(true),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(100, 200),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(100, 200)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(200, 300),
          aligned: true,
          new_block_id: BlockId(2),
          is_anchor: Some(false),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(300, 310),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(300, 310)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(310, 500),
          aligned: true,
          new_block_id: BlockId(3),
          is_anchor: Some(true),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(500, 600),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(500, 600)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(600, 900),
          aligned: true,
          new_block_id: BlockId(4),
          is_anchor: Some(false),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(900, 1000),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(900, 1000)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
      ]
    );
  }

  #[test]
  fn test_refine_intervals() -> Result<(), Report> {
    let (hits, block_length, bid) = example();
    let thr_len = 50;
    let intervals = extract_intervals(&hits, block_length, thr_len)?;

    assert_eq!(
      intervals,
      vec![
        PangraphInterval {
          interval: Interval::new(0, 100),
          aligned: true,
          new_block_id: BlockId(1),
          is_anchor: Some(true),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(100, 200),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(100, 200)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(200, 300),
          aligned: true,
          new_block_id: BlockId(2),
          is_anchor: Some(false),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(300, 500),
          aligned: true,
          new_block_id: BlockId(3),
          is_anchor: Some(true),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(500, 600),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(500, 600)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(600, 900),
          aligned: true,
          new_block_id: BlockId(4),
          is_anchor: Some(false),
          orientation: Some(Forward),
          cigar: None
        },
        PangraphInterval {
          interval: Interval::new(900, 1000),
          aligned: false,
          new_block_id: calculate_hash(bid, &Interval::new(900, 1000)),
          is_anchor: None,
          orientation: None,
          cigar: None
        },
      ]
    );

    Ok(())
  }
}
