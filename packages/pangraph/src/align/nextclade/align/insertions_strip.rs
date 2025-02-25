use crate::align::nextclade::alphabet::aa::Aa;
use crate::align::nextclade::alphabet::letter::{Letter, serde_deserialize_seq, serde_serialize_seq};
use crate::align::nextclade::alphabet::nuc::Nuc;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Insertion<T: Letter<T>> {
  pub pos: i32,

  #[serde(serialize_with = "serde_serialize_seq")]
  #[serde(deserialize_with = "serde_deserialize_seq")]
  pub ins: Vec<T>,
}

impl<T: Letter<T>> Insertion<T> {
  pub fn len(&self) -> usize {
    self.ins.len()
  }

  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }
}

pub type NucIns = Insertion<Nuc>;

/// Order nuc insertions by position, then length
impl<T: Letter<T>> Ord for Insertion<T> {
  fn cmp(&self, other: &Self) -> Ordering {
    (self.pos, self.ins.len()).cmp(&(other.pos, other.ins.len()))
  }
}

impl<T: Letter<T>> PartialOrd for Insertion<T> {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StripInsertionsResult<T: Letter<T>> {
  pub qry_seq: Vec<T>,
  pub insertions: Vec<Insertion<T>>,
}

pub fn insertions_strip<T: Letter<T>>(qry_seq: &[T], ref_seq: &[T]) -> StripInsertionsResult<T> {
  debug_assert_eq!(ref_seq.len(), qry_seq.len());

  let mut qry_stripped = Vec::<T>::with_capacity(ref_seq.len());
  let mut insertions = Vec::<Insertion<T>>::new();

  let mut insertion_start: i32 = -1;
  let mut ref_pos: i32 = -1;
  let mut current_insertion = Vec::<T>::with_capacity(16);

  for i in 0..ref_seq.len() {
    let c = ref_seq[i];

    if c.is_gap() {
      if current_insertion.is_empty() {
        // NOTE: by convention we set position of insertion to be the index of a character that precedes the insertion,
        // i.e. a position of reference nucleotide *after* which the insertion have happened.
        insertion_start = ref_pos;
      }
      current_insertion.push(qry_seq[i]);
    } else {
      qry_stripped.push(qry_seq[i]);
      ref_pos += 1;
      if !current_insertion.is_empty() {
        insertions.push(Insertion {
          pos: insertion_start,
          ins: current_insertion.clone(),
        });
        current_insertion.clear();
        insertion_start = -1;
      }
    }
  }

  if !current_insertion.is_empty() {
    insertions.push(Insertion {
      pos: insertion_start,
      ins: current_insertion.clone(),
    });
  }

  qry_stripped.shrink_to_fit();
  insertions.shrink_to_fit();
  insertions.sort();

  StripInsertionsResult {
    qry_seq: qry_stripped,
    insertions,
  }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AaIns {
  pub cds: String,
  pub pos: i32,

  #[serde(serialize_with = "serde_serialize_seq")]
  #[serde(deserialize_with = "serde_deserialize_seq")]
  pub ins: Vec<Aa>,
}

/// Order amino acid insertions by gene, position, then length
impl Ord for AaIns {
  fn cmp(&self, other: &Self) -> Ordering {
    (&self.cds, self.pos, self.ins.len()).cmp(&(&other.cds, other.pos, other.ins.len()))
  }
}

impl PartialOrd for AaIns {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}
