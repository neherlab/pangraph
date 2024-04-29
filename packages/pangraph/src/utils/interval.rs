use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct Interval {
  pub start: usize,
  pub end: usize,
}

impl Interval {
  pub fn new(start: usize, end: usize) -> Self {
    Self { start, end }
  }

  pub fn len(&self) -> usize {
    self.end.saturating_sub(self.start)
  }

  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }

  pub fn has_overlap_with(&self, other: &Interval) -> bool {
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
    assert!(!Interval::new(100, 200).has_overlap_with(&Interval::new(210, 390)));
    assert!(Interval::new(100, 220).has_overlap_with(&Interval::new(210, 390)));
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
