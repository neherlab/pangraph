use crate::{make_error, make_internal_report};
use eyre::Report;

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

/// Insert slice into vec at an index
/// Taken from: https://internals.rust-lang.org/t/add-vec-insert-slice-at-to-insert-the-content-of-a-slice-at-an-arbitrary-index/11008
pub fn insert_at_inplace<T: Clone>(vec: &mut Vec<T>, index: usize, slice: &[T]) {
  let len = vec.len();
  assert!(
    index <= len,
    "Attempted to insert outside of array boundaries: array size is {len}, index is {index}"
  );
  vec.reserve(slice.len());
  let mut v = vec.split_off(index);
  vec.extend_from_slice(slice);
  vec.append(&mut v);
}

#[cfg(test)]
mod tests {
  use super::*;
  
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_insert_at_inplace_general_case() {
    let mut vec = vec![1; 6];
    let slice = [0; 2];
    let index = 2;
    insert_at_inplace(&mut vec, index, &slice);
    assert_eq!(vec![1, 1, 0, 0, 1, 1, 1, 1], vec);
  }

  #[rstest]
  fn test_insert_at_inplace_append() {
    let mut vec = vec![1; 6];
    let slice = [0; 2];
    let index = 6;
    insert_at_inplace(&mut vec, index, &slice);
    assert_eq!(vec![1, 1, 1, 1, 1, 1, 0, 0], vec);
  }

  #[rstest]
  fn test_insert_at_inplace_prepend() {
    let mut vec = vec![1; 6];
    let slice = [0; 2];
    let index = 0;
    insert_at_inplace(&mut vec, index, &slice);
    assert_eq!(vec![0, 0, 1, 1, 1, 1, 1, 1], vec);
  }
}
