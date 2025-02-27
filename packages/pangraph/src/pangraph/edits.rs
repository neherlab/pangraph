use crate::io::seq::{complement, reverse_complement};
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use crate::utils::interval::Interval;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::ops::Range;

#[cfg(any(test, debug_assertions))]
use eyre::eyre;

#[must_use]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, JsonSchema)]
pub struct Sub {
  pub pos: usize,
  pub alt: AsciiChar,
}

impl Sub {
  pub fn new(pos: usize, alt: impl Into<AsciiChar>) -> Self {
    Self { pos, alt: alt.into() }
  }

  pub fn reverse_complement(&self, len: usize) -> Result<Self, Report> {
    Ok(Self {
      pos: len - self.pos - 1,
      alt: complement(&self.alt)?,
    })
  }

  pub fn shift(&self, shift: isize) -> Self {
    Self {
      pos: (self.pos as isize + shift) as usize,
      alt: self.alt,
    }
  }
}

#[must_use]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, JsonSchema)]
pub struct Del {
  pub pos: usize,
  pub len: usize,
}

impl Del {
  pub fn new(pos: usize, len: usize) -> Self {
    Self { pos, len }
  }

  pub fn end(&self) -> usize {
    self.pos + self.len
  }

  pub fn range(&self) -> Range<usize> {
    self.pos..self.end()
  }

  pub fn interval(&self) -> Interval {
    Interval::new(self.pos, self.end())
  }

  pub fn reverse_complement(&self, len: usize) -> Self {
    Self {
      pos: len - self.pos - self.len,
      len: self.len,
    }
  }

  pub fn shift(&self, shift: isize) -> Self {
    Self {
      pos: (self.pos as isize + shift) as usize,
      len: self.len,
    }
  }

  pub fn contains(&self, pos: usize) -> bool {
    self.range().contains(&pos)
  }
}

#[must_use]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Ord, PartialOrd, Hash, JsonSchema)]
pub struct Ins {
  pub pos: usize,
  pub seq: Seq,
}

impl Ins {
  pub fn new(pos: usize, seq: impl Into<Seq>) -> Self {
    Self { pos, seq: seq.into() }
  }

  pub fn reverse_complement(&self, len: usize) -> Result<Self, Report> {
    Ok(Self {
      pos: len - self.pos,
      seq: reverse_complement(&self.seq)?,
    })
  }

  pub fn shift(&self, shift: isize) -> Self {
    Self {
      pos: self.pos.saturating_add_signed(shift),
      seq: self.seq.clone(),
    }
  }
}

#[must_use]
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq, Hash, JsonSchema)]
pub struct Edit {
  pub subs: Vec<Sub>,
  pub dels: Vec<Del>,
  pub inss: Vec<Ins>,
}

impl Edit {
  pub fn empty() -> Self {
    Self::default()
  }

  pub fn is_empty(&self) -> bool {
    self.subs.is_empty() && self.dels.is_empty() && self.inss.is_empty()
  }

  /// Construct edit which consists of a deletion of length `len`
  pub fn deleted(len: usize) -> Self {
    Self {
      subs: vec![],
      dels: vec![Del::new(0, len)],
      inss: vec![],
    }
  }

  pub fn new(inss: impl Into<Vec<Ins>>, dels: impl Into<Vec<Del>>, subs: impl Into<Vec<Sub>>) -> Self {
    Self {
      subs: subs.into(),
      dels: dels.into(),
      inss: inss.into(),
    }
  }

  pub fn reverse_complement(&self, len: usize) -> Result<Self, Report> {
    let subs = self
      .subs
      .iter()
      .map(|s| s.reverse_complement(len))
      .collect::<Result<Vec<_>, Report>>()?;

    let dels = self.dels.iter().map(|d| d.reverse_complement(len)).collect_vec();

    let inss = self
      .inss
      .iter()
      .map(|i| i.reverse_complement(len))
      .collect::<Result<Vec<_>, Report>>()?;

    Ok(Self { subs, dels, inss })
  }

  pub fn shift(&self, shift: isize) -> Self {
    Self {
      inss: self.inss.iter().map(|i| i.shift(shift)).collect(),
      dels: self.dels.iter().map(|d| d.shift(shift)).collect(),
      subs: self.subs.iter().map(|s| s.shift(shift)).collect(),
    }
  }

  pub fn concat(&self, next: &Self) -> Self {
    let mut inss = self.inss.clone();
    let mut dels = self.dels.clone();
    let mut subs = self.subs.clone();

    // if two insertions have the same position,
    // concatenate self + next
    for ins in &next.inss {
      if let Some(prev) = inss.iter_mut().find(|i| i.pos == ins.pos) {
        prev.seq.extend_seq(&ins.seq);
      } else {
        inss.push(ins.clone());
      }
    }

    dels.extend(next.dels.iter().cloned());
    subs.extend(next.subs.iter().cloned());
    Self { subs, dels, inss }
  }

  /// Apply the edits to the reference to obtain the query sequence
  pub fn apply(&self, reff: impl Into<Seq>) -> Result<Seq, Report> {
    let mut qry = reff.into();

    for sub in &self.subs {
      qry[sub.pos] = sub.alt;
    }

    for del in &self.dels {
      // Replace deleted nucs with character `-`, to avoid frame shift
      for pos in del.range() {
        qry[pos] = AsciiChar(b'-');
      }
    }

    for ins in self.inss.iter().sorted().rev() {
      qry.insert_seq(ins.pos, &ins.seq);
    }

    // Strip gaps introduced when applying deletions
    qry.retain(|c| c != &AsciiChar(b'-'));

    Ok(qry)
  }

  pub fn apply_aligned(&self, reff: impl Into<Seq>) -> Result<Seq, Report> {
    let mut qry = reff.into();

    for sub in &self.subs {
      qry[sub.pos] = sub.alt;
    }

    for del in &self.dels {
      for pos in del.range() {
        qry[pos] = AsciiChar(b'-');
      }
    }

    Ok(qry)
  }

  /// Check if the alignment is empty, i.e. after applying the edits to the
  /// consensus sequence, the resulting sequence is empty.
  pub fn is_empty_alignment(&self, consensus: &Seq) -> bool {
    let cons_len = consensus.len();
    // if there are insertions, the alignment is not empty
    let insertions_len = self.inss.iter().map(|i| i.seq.len()).sum::<usize>();
    if insertions_len > 0 {
      return false;
    }
    // if the total length of deletions is less than the length of the consensus
    // sequence, the alignment is not empty
    let deletions_len = self.dels.iter().map(|d| d.len).sum::<usize>();
    if deletions_len < cons_len {
      return false;
    }
    // otherwise, apply the edits and check if the resulting sequence is empty
    let seq_len = self.apply(consensus).unwrap().len();
    seq_len == 0
  }

  /// cumulative length of leading deletions, if present
  /// leading deletions are deletions that start at the beginning of the sequence
  /// Nb: assumes that deletions are not adjacent
  pub fn leading_deletions(&self) -> usize {
    self.dels.iter().filter(|d| d.pos == 0).map(|d| d.len).sum()
  }

  /// cumulative length of trailing deletions, if present
  pub fn trailing_deletions(&self, consensus_len: usize) -> usize {
    self
      .dels
      .iter()
      .filter(|d| d.end() == consensus_len)
      .map(|d| d.len)
      .sum()
  }

  /// total number of internal deletions
  pub fn internal_deletions(&self, consensus_len: usize) -> usize {
    let leading = self.leading_deletions();
    let trailing = self.trailing_deletions(consensus_len);
    let total = self.dels.iter().map(|d| d.len).sum::<usize>();
    total - leading - trailing
  }

  /// cumulative length of leading insertions, if present
  pub fn leading_insertions(&self) -> usize {
    self.inss.iter().filter(|i| i.pos == 0).map(|i| i.seq.len()).sum()
  }

  /// cumulative length of trailing insertions, if present
  pub fn trailing_insertions(&self, consensus_len: usize) -> usize {
    self
      .inss
      .iter()
      .filter(|i| i.pos == consensus_len)
      .map(|i| i.seq.len())
      .sum()
  }

  /// total number of internal insertions
  pub fn internal_insertions(&self, consensus_len: usize) -> usize {
    let leading = self.leading_insertions();
    let trailing = self.trailing_insertions(consensus_len);
    let total = self.inss.iter().map(|i| i.seq.len()).sum::<usize>();
    total - leading - trailing
  }

  /// Returns the number of aligned positions after position `p` in the consensus sequence.
  pub fn aligned_count_after(&self, p: usize, cons_len: usize) -> usize {
    let total_positions = cons_len.saturating_sub(p);
    let deletion_overlap: usize = self
      .dels
      .iter()
      .filter_map(|d| {
        // If the deletion is completely before p ignore it.
        if d.end() <= p {
          None
        } else {
          // The overlapping region is from max(p, d.pos) to min(cons_len, d.end())
          let start = p.max(d.pos);
          Some(d.end() - start)
        }
      })
      .sum();
    total_positions.saturating_sub(deletion_overlap)
  }

  /// Given an alignment, returns the average displacement of the sequence compared
  /// to the consensus. This is the average absolute difference between the position
  /// of the aligned nucleotide on the query (after applying edits) and its position
  /// on the consensus.
  ///
  /// ref : xxxxxxxxxxxxxx---
  /// qry : ---xxxxxxxxxxxxxx
  /// average displacement = +3
  ///
  /// ref : ---xxxxxxxxxxxxx
  /// qry : xxxxxxxxxxxxx---
  /// average displacement = -3
  pub fn average_displacement(&self, cons_len: usize) -> Option<i32> {
    // Total number of aligned positions.
    let aligned_count = self.aligned_count_after(0, cons_len);
    if aligned_count == 0 {
      return None;
    }

    let mut total_shift = 0_i32;
    // Insertions push the query to the right, so they cause a negative displacement.
    for ins in &self.inss {
      let p = ins.pos;
      let count_after = self.aligned_count_after(p, cons_len);
      total_shift -= (ins.seq.len() * count_after) as i32;
    }
    // Deletions pull the query to the left (i.e. the query lacks bases),
    // so they contribute a positive displacement.
    for d in &self.dels {
      let p = d.pos;
      let count_after = self.aligned_count_after(p, cons_len);
      total_shift += (d.len * count_after) as i32;
    }
    Some((total_shift as f64 / aligned_count as f64).round() as i32)
  }

  #[cfg(any(test, debug_assertions))]
  pub fn sanity_check(&self, len: usize) -> Result<(), Report> {
    let block_interval = Interval::new(0, len);
    // === substitution checks ===
    for sub in &self.subs {
      if !block_interval.contains(sub.pos) {
        return Err(eyre!(
          "Substitution position {} is out of bounds for sequence of length {}",
          sub.pos,
          len
        ));
      }
      if sub.alt == AsciiChar(b'-') {
        return Err(eyre!("Substitution with deletion character '-' is not allowed"));
      }
    }

    // check that there are no two substitutions at the same position
    for i in 0..self.subs.len() {
      for j in i + 1..self.subs.len() {
        if self.subs[i].pos == self.subs[j].pos {
          return Err(eyre!(
            "Substitution {:?} overlaps with substitution {:?}",
            self.subs[i],
            self.subs[j]
          ));
        }
      }
    }

    // check that substitutions do not overlap with deletions
    for sub in &self.subs {
      for del in &self.dels {
        if del.contains(sub.pos) {
          return Err(eyre!("Substitution {:?} overlaps with deletion {:?}", sub, del));
        }
      }
    }

    // === deletion checks ===
    for del in &self.dels {
      if del.len == 0 {
        return Err(eyre!("Deletion {:?} has length 0", del));
      }

      if !block_interval.contains(del.pos) {
        return Err(eyre!(
          "Deletion {:?} is out of bounds for sequence of length {}",
          del,
          len
        ));
      }

      if del.end() > len {
        return Err(eyre!(
          "Deletion {:?} is out of bounds for sequence of length {}",
          del,
          len
        ));
      }
    }

    // check that deletions are non-overlapping
    // for all pairs of deletions
    for i in 0..self.dels.len() {
      for j in i + 1..self.dels.len() {
        let del_i = &self.dels[i];
        let del_j = &self.dels[j];
        if del_i.interval().has_overlap_with(&del_j.interval()) {
          return Err(eyre!("Deletion {:?} overlaps with deletion {:?}", del_i, del_j));
        }
      }
    }

    // === insertion checks ===
    for ins in &self.inss {
      if ins.pos > len {
        return Err(eyre!(
          "Insertion position {} is out of bounds for sequence of length {}",
          ins.pos,
          len
        ));
      }
    }

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_edits_apply_simple_case() {
    //                 1
    //       01234567890
    let r = "ACCTGGCTTT";
    let q = "AGCCTGGTAT";

    let inss = vec![Ins::new(1, "G")];
    let dels = vec![Del::new(6, 1)];
    let subs = vec![Sub::new(8, 'A')];
    let edits = Edit { subs, dels, inss };

    let actual = edits.apply(r).unwrap();
    assert_eq!(q, &actual);
  }

  #[rstest]
  fn test_edits_apply_sub() {
    //       0123456789
    let r = "ACCTGGCTTT";
    let q = "ACCAGGCTTT";
    //          s

    let edits = Edit {
      subs: vec![Sub::new(3, 'A')],
      dels: vec![],
      inss: vec![],
    };

    let actual = edits.apply(r).unwrap();
    assert_eq!(q, &actual);
  }

  #[rstest]
  fn test_edits_apply_del() {
    //       0123456789
    let r = "ACCTGGCTTT";
    let q = "ACCGCTTT";
    //          d

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(3, 2)],
      inss: vec![],
    };

    let actual = edits.apply(r).unwrap();
    assert_eq!(q, &actual);
  }

  #[rstest]
  fn test_edits_apply_ins() {
    //
    //       0123456789
    let r = "ACCTGGCTTT";
    let q = "ACCTACGGCTTT";
    //           ii

    let edits = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(4, "AC")],
    };

    let actual = edits.apply(r).unwrap();
    assert_eq!(q, &actual);
  }

  #[rstest]
  fn test_empty_alignment() {
    let consensus = "ACGT";
    let edits = Edit::empty();
    assert!(!edits.is_empty_alignment(&Seq::from(consensus)));

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 4)],
      inss: vec![Ins::new(1, "A")],
    };
    assert!(!edits.is_empty_alignment(&Seq::from(consensus)));

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 4)],
      inss: vec![],
    };
    assert!(edits.is_empty_alignment(&Seq::from(consensus)));
  }

  #[rstest]
  fn test_leading_deletions() {
    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 3), Del::new(6, 2)],
      inss: vec![],
    };
    assert_eq!(edits.leading_deletions(), 3);

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(1, 3)],
      inss: vec![],
    };
    assert_eq!(edits.leading_deletions(), 0);
  }

  #[rstest]
  fn test_trailing_deletions() {
    let consensus_len = 10;
    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(8, 2), Del::new(0, 3)],
      inss: vec![],
    };
    assert_eq!(edits.trailing_deletions(consensus_len), 2);

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(4, 3)],
      inss: vec![],
    };
    assert_eq!(edits.trailing_deletions(consensus_len), 0);
  }

  #[rstest]
  fn test_internal_deletions() {
    let consensus_len = 10;
    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 2), Del::new(4, 2), Del::new(8, 2)],
      inss: vec![],
    };
    assert_eq!(edits.internal_deletions(consensus_len), 2);

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(1, 3), Del::new(5, 2)],
      inss: vec![],
    };
    assert_eq!(edits.internal_deletions(consensus_len), 5);
  }

  #[rstest]
  fn test_leading_insertions() {
    let edits = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(0, "AAA"), Ins::new(5, "GGG")],
    };
    assert_eq!(edits.leading_insertions(), 3);

    let edits = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(1, "AAA")],
    };
    assert_eq!(edits.leading_insertions(), 0);
  }

  #[rstest]
  fn test_trailing_insertions() {
    let consensus_len = 10;
    let edits = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(10, "TTT"), Ins::new(0, "AAAA")],
    };
    assert_eq!(edits.trailing_insertions(consensus_len), 3);

    let edits = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(5, "TTT")],
    };
    assert_eq!(edits.trailing_insertions(consensus_len), 0);
  }

  #[rstest]
  fn test_internal_insertions() {
    let consensus_len = 10;
    let edits = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(0, "AAA"), Ins::new(5, "GGG"), Ins::new(10, "TTT")],
    };
    assert_eq!(edits.internal_insertions(consensus_len), 3);

    let edits = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(1, "AAA"), Ins::new(5, "GGG")],
    };
    assert_eq!(edits.internal_insertions(consensus_len), 6);
  }

  #[rstest]
  fn test_aligned_count_no_deletions() {
    // When there are no deletions, the aligned count is just cons_len minus p.
    let edit = Edit::empty();
    let cons_len = 10;
    assert_eq!(edit.aligned_count_after(0, cons_len), 10);
    assert_eq!(edit.aligned_count_after(5, cons_len), 5);
    assert_eq!(edit.aligned_count_after(10, cons_len), 0);
  }

  #[rstest]
  fn test_aligned_count_single_deletion() {
    // A deletion covering positions 3 and 4.
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(3, 2)],
      inss: vec![],
    };
    let cons_len = 10;
    assert_eq!(edit.aligned_count_after(0, cons_len), 8);
    assert_eq!(edit.aligned_count_after(2, cons_len), 6);
    assert_eq!(edit.aligned_count_after(4, cons_len), 5);
    assert_eq!(edit.aligned_count_after(5, cons_len), 5);
    assert_eq!(edit.aligned_count_after(10, cons_len), 0);
  }

  #[rstest]
  fn test_aligned_count_after_multiple_deletions() {
    // Two deletions: one covering positions 3..7 and one covering positions 10..13.
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(3, 4), Del::new(10, 3)],
      inss: vec![],
    };
    let cons_len = 20;
    assert_eq!(edit.aligned_count_after(0, cons_len), 13);
    assert_eq!(edit.aligned_count_after(5, cons_len), 10);
    assert_eq!(edit.aligned_count_after(12, cons_len), 7);
    assert_eq!(edit.aligned_count_after(13, cons_len), 7);
    assert_eq!(edit.aligned_count_after(17, cons_len), 3);
  }

  #[rstest]
  fn test_average_displacement_no_edits() {
    // No edits: expect zero displacement.
    let edit = Edit::empty();
    let cons_len = 10;
    assert_eq!(edit.average_displacement(cons_len), Some(0));
  }

  #[rstest]
  fn test_average_displacement_insertion() {
    // Insertion at the beginning: all positions receive a shift of -2.
    let edit = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(0, "AA")],
    };
    let cons_len = 10;
    assert_eq!(edit.average_displacement(cons_len), Some(-2));

    // Insertion at the end: all positions receive a shift of 0.
    let edit = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(10, "AA")],
    };
    let cons_len = 10;
    assert_eq!(edit.average_displacement(cons_len), Some(0));
  }

  #[rstest]
  fn test_average_displacement_deletion() {
    // A leading deletion of length 2: all positions receive a shift of 2.
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(2, 2)],
      inss: vec![],
    };
    let cons_len = 10;
    assert_eq!(edit.average_displacement(cons_len), Some(2));

    // A trailing deletion of length 2: all positions receive a shift of 0.
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(8, 2)],
      inss: vec![],
    };
    let cons_len = 10;
    assert_eq!(edit.average_displacement(cons_len), Some(0));
  }

  #[rstest]
  fn test_average_displacement_ins_and_del() {
    let cons_len = 10;

    // A leading deletion of length 3 and insertion of length 2: all positions receive a shift of 1.
    //     012  3456789
    // R = xxx--xxxxxxx
    // Q = ---xxxxxxxxx
    //        012345678
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 3)],
      inss: vec![Ins::new(3, "AA")],
    };
    assert_eq!(edit.average_displacement(cons_len), Some(1));

    // An internal insertion of length 4
    //     01234    56789
    // R = xxxxx----xxxxx
    // Q = xxxxxxxxxxxxxx
    //     01234567890123
    let edit = Edit {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(4, "AAAA")],
    };
    assert_eq!(edit.average_displacement(cons_len), Some(-2));

    // An internal deletion of length 3
    //     0123456789
    // R = xxxxxxxxxx
    // Q = xxxx---xxx
    //     0123   456
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(4, 3)],
      inss: vec![],
    };
    assert_eq!(edit.average_displacement(cons_len), Some(1));

    // Two deletions and insertions
    //        0123  45678901
    // R = ---xxxx--xxxxxxxx
    // Q = xxxxx--xxxx---xxx
    //     01234  5678   901
    let cons_len = 12;
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(2, 2), Del::new(6, 3)],
      inss: vec![Ins::new(0, "AAA"), Ins::new(4, "AA")],
    };
    assert_eq!(edit.average_displacement(cons_len), Some(-2));
  }

  #[rstest]
  fn test_average_displacement_full_deletion() {
    // A deletion that removes the entire consensus results in no aligned positions,
    // so the average displacement is defined as 0.
    let edit = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 10)],
      inss: vec![],
    };
    let cons_len = 10;
    assert_eq!(edit.average_displacement(cons_len), None);
  }
}
