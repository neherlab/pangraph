#[derive(Debug)]
pub struct Interval {
  start: i32,
  end: i32,
}

impl Interval {
  fn overlaps(&self, i: &Interval) -> bool {
    // self |-----|       or   |--------|
    // i        |-----|          |----|
    if (i.start >= self.start) && (i.start <= self.end) {
      return true;
    }
    // self     |-----|
    // i     |-----|
    else if (i.end >= self.start) && (i.end <= self.end) {
      return true;
    }
    // self     |---|
    // i     |----------|
    else if (self.start >= i.start) && (self.start <= i.end) {
      return true;
    }
    false
  }
}

fn no_overlap(intervals: &[Interval], candidate: &Interval) -> bool {
  !intervals.iter().any(|interval| interval.overlaps(candidate))
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_overlap() {
    assert!(!Interval { start: 100, end: 200 }.overlaps(&Interval { start: 210, end: 390 }));
    assert!(Interval { start: 100, end: 220 }.overlaps(&Interval { start: 210, end: 390 }));
  }

  #[test]
  fn test_overlap_no_overlap() {
    assert!(no_overlap(
      &[Interval { start: 100, end: 200 }, Interval { start: 300, end: 400 }],
      &Interval { start: 210, end: 290 }
    ));
    assert!(!no_overlap(
      &[Interval { start: 100, end: 200 }, Interval { start: 300, end: 400 }],
      &Interval { start: 210, end: 390 }
    ));
  }
}
