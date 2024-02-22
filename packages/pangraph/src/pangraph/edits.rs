use crate::utils::collections::insert_at_inplace;
use eyre::Report;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::ops::Range;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Sub {
  pos: usize,
  alt: u8,
}

impl Sub {
  pub fn new(pos: usize, alt: u8) -> Self {
    Self { pos, alt }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Del {
  pos: usize,
  len: usize,
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
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Ins {
  pos: usize,
  seq: String,
}

impl Ins {
  pub fn new(pos: usize, seq: impl AsRef<str>) -> Self {
    Self {
      pos,
      seq: seq.as_ref().to_owned(),
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Edits {
  subs: Vec<Sub>,
  dels: Vec<Del>,
  inss: Vec<Ins>,
}

impl Edits {
  /// Apply the edits to the reference to obtain the query sequence
  pub fn apply(&self, reff: impl AsRef<str>) -> Result<String, Report> {
    let mut qry: Vec<u8> = reff.as_ref().to_owned().into();

    for sub in &self.subs {
      qry[sub.pos] = sub.alt;
    }

    for del in &self.dels {
      qry.drain(del.range());
    }

    for ins in &self.inss {
      insert_at_inplace(&mut qry, ins.pos, ins.seq.as_bytes());
    }

    let qry = String::from_utf8(qry)?;
    Ok(qry)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
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
    let subs = vec![Sub::new(8, b'A')];
    let edits = Edits { subs, dels, inss };

    let actual = edits.apply(r).unwrap();
    assert_eq!(q, actual);
  }

  #[rstest]
  fn test_edits_apply_sub() {
    //       0123456789
    let r = "ACCTGGCTTT";
    let q = "ACCAGGCTTT";
    //          s

    let edits = Edits {
      subs: vec![Sub::new(3, b'A')],
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

    let edits = Edits {
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

    let edits = Edits {
      subs: vec![],
      dels: vec![],
      inss: vec![Ins::new(4, "AC")],
    };

    let actual = edits.apply(r).unwrap();
    assert_eq!(q, actual);
  }
}
