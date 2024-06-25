use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct Interval {
  pub start: usize,
  pub end: usize,
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
}

pub fn have_no_overlap(intervals: &[Interval], candidate: &Interval) -> bool {
  !intervals.iter().any(|interval| interval.has_overlap_with(candidate))
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
}
