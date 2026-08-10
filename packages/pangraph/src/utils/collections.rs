use crate::{make_error, make_internal_report};
use eyre::Report;
use std::collections::{BTreeMap, HashSet};
use std::hash::Hash;

pub fn concat_to_vec<T: Clone>(x: &[T], y: &[T]) -> Vec<T> {
  [x, y].into_iter().flatten().cloned().collect()
}

pub fn first<T>(arr: &[T]) -> Result<&T, Report> {
  arr
    .first()
    .ok_or_else(|| make_internal_report!("When attempted to retrieve the first element: Array is empty"))
}

pub fn last<T>(arr: &[T]) -> Result<&T, Report> {
  arr
    .last()
    .ok_or_else(|| make_internal_report!("When attempted to retrieve the last element: Array is empty"))
}

// Ensure the vector has exactly one element and return it
pub fn remove_exactly_one<T>(mut elems: Vec<T>) -> Result<T, Report> {
  match elems.len() {
    0 => make_error!("Expected exactly one element, but found none"),
    1 => Ok(elems.remove(0)),
    _ => make_error!("Expected exactly one element, but found: {}", elems.len()),
  }
}

pub fn has_duplicates<T: Eq + Hash, I: IntoIterator<Item = T>>(iter: I) -> bool {
  let mut seen = HashSet::new();
  iter.into_iter().any(|item| !seen.insert(item))
}

/// Returns the values that occur more than once in the given iterator, in sorted order.
/// Each duplicated value is reported once, no matter how many times it occurs.
pub fn find_duplicates<T: Ord, I: IntoIterator<Item = T>>(iter: I) -> Vec<T> {
  let mut counts: BTreeMap<T, usize> = BTreeMap::new();
  for item in iter {
    *counts.entry(item).or_insert(0) += 1;
  }
  counts
    .into_iter()
    .filter(|(_, count)| *count > 1)
    .map(|(item, _)| item)
    .collect()
}
