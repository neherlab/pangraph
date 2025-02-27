use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter};
use std::ops::Range;

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct Interval {
  pub start: usize,
  pub end: usize,
}

impl Display for Interval {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "[{}, {})", self.start, self.end)
  }
}

impl Interval {
  pub fn new(start: usize, end: usize) -> Self {
    let this = Self { start, end };
    debug_assert!(Interval::is_valid(&this));
    this
  }

  pub fn is_valid(&self) -> bool {
    self.end >= self.start
  }

  pub fn len(&self) -> usize {
    debug_assert!(self.is_valid());
    self.end.saturating_sub(self.start)
  }

  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }

  pub fn contains(&self, pos: usize) -> bool {
    self.start <= pos && pos < self.end
  }

  pub fn has_overlap_with(&self, other: &Interval) -> bool {
    debug_assert!(self.is_valid());
    debug_assert!(other.is_valid());
    self.end > other.start && self.start < other.end
  }

  pub fn to_range(&self) -> Range<usize> {
    Range {
      start: self.start,
      end: self.end,
    }
  }
}

pub fn have_no_overlap(intervals: &[Interval], candidate: &Interval) -> bool {
  !intervals.iter().any(|interval| interval.has_overlap_with(candidate))
}

pub fn positions_to_intervals(positions: &[usize]) -> Vec<Interval> {
  if positions.is_empty() {
    return Vec::new();
  }
  let mut sorted = positions.to_vec();
  sorted.sort_unstable();

  let mut intervals = Vec::new();
  let mut start = sorted[0];
  let mut end = sorted[0];

  for &pos in sorted.iter().skip(1) {
    // If the position is contiguous (or a duplicate), extend the current interval.
    if pos > end + 1 {
      intervals.push(Interval::new(start, end + 1));
      start = pos;
    }
    end = pos;
  }
  intervals.push(Interval::new(start, end + 1));
  intervals
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_overlap() {
    // No overlap: self is completely before other
    // [1, 5) [5, 10)
    assert!(!Interval::new(1, 5).has_overlap_with(&Interval::new(5, 10)));
    assert!(!Interval::new(5, 10).has_overlap_with(&Interval::new(1, 5)));

    // No overlap: self is completely after other
    // [10, 15) [5, 10)
    assert!(!Interval::new(10, 15).has_overlap_with(&Interval::new(5, 10)));
    assert!(!Interval::new(5, 10).has_overlap_with(&Interval::new(10, 15)));

    // Partial overlap: self starts before and ends within other
    // [1, 7) [5, 10)
    assert!(Interval::new(1, 7).has_overlap_with(&Interval::new(5, 10)));
    assert!(Interval::new(5, 10).has_overlap_with(&Interval::new(1, 7)));

    // Partial overlap: self starts within and ends after other
    // [5, 12) [5, 10)
    assert!(Interval::new(5, 12).has_overlap_with(&Interval::new(5, 10)));
    assert!(Interval::new(5, 10).has_overlap_with(&Interval::new(5, 12)));

    // Full overlap: self completely covers other
    // [1, 15) [5, 10)
    assert!(Interval::new(1, 15).has_overlap_with(&Interval::new(5, 10)));
    assert!(Interval::new(5, 10).has_overlap_with(&Interval::new(1, 15)));

    // Full overlap: other completely covers self
    // [5, 10) [1, 15)
    assert!(Interval::new(5, 10).has_overlap_with(&Interval::new(1, 15)));
    assert!(Interval::new(1, 15).has_overlap_with(&Interval::new(5, 10)));

    // Overlap at the boundary where self.end == other.start
    // [1, 5) [5, 10)
    assert!(!Interval::new(1, 5).has_overlap_with(&Interval::new(5, 10)));
    assert!(!Interval::new(5, 10).has_overlap_with(&Interval::new(1, 5)));

    // Overlap at the boundary where self.start == other.end
    // [5, 10) [1, 5)
    assert!(!Interval::new(5, 10).has_overlap_with(&Interval::new(1, 5)));
    assert!(!Interval::new(1, 5).has_overlap_with(&Interval::new(5, 10)));
  }

  #[test]
  fn test_no_overlap() {
    assert!(have_no_overlap(
      &[Interval::new(100, 200), Interval::new(300, 400)],
      &Interval::new(210, 290)
    ));
    assert!(!have_no_overlap(
      &[Interval::new(100, 200), Interval::new(300, 400)],
      &Interval::new(210, 390)
    ));
  }

  #[test]
  fn test_from_position_list_empty() {
    let positions: Vec<usize> = vec![];
    let intervals = positions_to_intervals(&positions);
    assert!(intervals.is_empty());
  }

  #[test]
  fn test_from_position_list_single() {
    let positions = vec![5];
    let intervals = positions_to_intervals(&positions);
    assert_eq!(intervals, vec![Interval::new(5, 6)]);
  }

  #[test]
  fn test_from_position_list_contiguous() {
    let positions = vec![1, 2, 3, 4, 5];
    let intervals = positions_to_intervals(&positions);
    // [1, 5] becomes [1, 6) as our intervals are half-open.
    assert_eq!(intervals, vec![Interval::new(1, 6)]);
  }

  #[test]
  fn test_from_position_list_non_contiguous() {
    let positions = vec![1, 3, 5];
    let intervals = positions_to_intervals(&positions);
    // Each position stands alone.
    assert_eq!(
      intervals,
      vec![Interval::new(1, 2), Interval::new(3, 4), Interval::new(5, 6)]
    );
  }

  #[test]
  fn test_from_position_list_unsorted() {
    let positions = vec![10, 1, 2, 3, 20, 21];
    let intervals = positions_to_intervals(&positions);
    // After sorting: [1,2,3,10,20,21]
    assert_eq!(
      intervals,
      vec![
        Interval::new(1, 4),   // from 1 to 3 -> [1,4)
        Interval::new(10, 11), // singleton [10,11)
        Interval::new(20, 22)  // 20 and 21 merge into [20,22)
      ]
    );
  }

  #[test]
  fn test_from_position_list_duplicates() {
    let positions = vec![5, 5, 5, 6, 7, 7, 8];
    let intervals = positions_to_intervals(&positions);
    // Duplicates should be merged into a single contiguous interval.
    assert_eq!(intervals, vec![Interval::new(5, 9)]);
  }
}
