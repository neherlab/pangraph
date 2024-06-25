use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::strand::Strand;
use crate::utils::interval::Interval;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PangraphInterval {
  pub interval: Interval,
  pub aligned: bool,
  pub new_block_id: BlockId,
  pub is_anchor: bool,
  pub orientation: Strand,
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

  pub fn has_overlap_with(&self, other: &PangraphInterval) -> bool {
    self.interval.has_overlap_with(&other.interval)
  }

  pub fn insertion_overlap(&self, ins_pos: usize, block_len: usize) -> bool {
    self.interval.contains(ins_pos) || (ins_pos == block_len && self.interval.end == block_len)
  }
}

pub fn have_no_overlap(intervals: &[PangraphInterval], candidate: &PangraphInterval) -> bool {
  !intervals.iter().any(|interval| interval.has_overlap_with(candidate))
}
