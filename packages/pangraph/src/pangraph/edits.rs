use crate::io::seq::{complement, reverse_complement};
use crate::utils::collections::insert_at_inplace;
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
  pub alt: char,
}

impl Sub {
  pub fn new(pos: usize, alt: char) -> Self {
    Self { pos, alt }
  }

  pub fn reverse_complement(&self, len: usize) -> Result<Self, Report> {
    Ok(Self {
      pos: len - self.pos - 1,
      alt: complement(self.alt)?,
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
  pub seq: String,
}

impl Ins {
  pub fn new(pos: usize, seq: impl AsRef<str>) -> Self {
    Self {
      pos,
      seq: seq.as_ref().to_owned(),
    }
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
        prev.seq.push_str(&ins.seq);
      } else {
        inss.push(ins.clone());
      }
    }

    dels.extend(next.dels.iter().cloned());
    subs.extend(next.subs.iter().cloned());
    Self { subs, dels, inss }
  }

  /// Apply the edits to the reference to obtain the query sequence
  pub fn apply(&self, reff: impl AsRef<str>) -> Result<String, Report> {
    // TODO: decide whether it's best to use chars, bytes of something else entirely
    let mut qry: Vec<char> = reff.as_ref().chars().collect_vec();

    for sub in &self.subs {
      qry[sub.pos] = sub.alt;
    }

    for del in &self.dels {
      // Replace deleted nucs with character `-`, to avoid frame shift
      for pos in del.range() {
        qry[pos] = '-';
      }
    }

    for ins in self.inss.iter().sorted().rev() {
      // TODO: avoid copy
      let seq = ins.seq.chars().collect_vec();
      insert_at_inplace(&mut qry, ins.pos, &seq);
    }

    // Strip gaps introduced when applying deletions
    qry.retain(|c| c != &'-');

    let qry = String::from_iter(&qry);
    Ok(qry)
  }

  pub fn apply_aligned(&self, reff: impl AsRef<str>) -> Result<String, Report> {
    // TODO: decide whether it's best to use chars, bytes of something else entirely
    let mut qry: Vec<char> = reff.as_ref().chars().collect_vec();

    for sub in &self.subs {
      qry[sub.pos] = sub.alt;
    }

    for del in &self.dels {
      for pos in del.range() {
        qry[pos] = '-';
      }
    }

    let qry = String::from_iter(&qry);
    Ok(qry)
  }

  /// Check if the alignment is empty, i.e. after applying the edits to the
  /// consensus sequence, the resulting sequence is empty.
  pub fn is_empty_alignment(&self, consensus: impl AsRef<str>) -> bool {
    let cons_len = consensus.as_ref().len();
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
    return seq_len == 0;
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
      if sub.alt == '-' {
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
    assert_eq!(q, actual);
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
    assert_eq!(q, actual);
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
    assert_eq!(q, actual);
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
    assert_eq!(q, actual);
  }

  #[rstest]
  fn test_empty_alignment() {
    let consensus = "ACGT";
    let edits = Edit::empty();
    assert!(!edits.is_empty_alignment(consensus));

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 4)],
      inss: vec![Ins::new(1, "A")],
    };
    assert!(!edits.is_empty_alignment(consensus));

    let edits = Edit {
      subs: vec![],
      dels: vec![Del::new(0, 4)],
      inss: vec![],
    };
    assert!(edits.is_empty_alignment(consensus));
  }
}
