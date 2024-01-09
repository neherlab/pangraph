use crate::distance::mash::minimizer::{minimizers_sketch, MinimizersParams};
use crate::pangraph::pangraph::{Pangraph, PangraphBlock, PangraphPath};
use itertools::Itertools;
use ndarray::{array, Array, Array2};
use serde::{Deserialize, Serialize};

/// Compute the pairwise distance between all input graphs.
/// Distance is the set distance between minimizers.
/// Linear-time algorithm using hash collisions.
pub fn mash_distance(graphs: &[Pangraph], params: &MinimizersParams) -> Array2<u64> {
  if graphs.is_empty() {
    return array![[]];
  }

  let sequences = graphs
    .iter()
    .flat_map(|graph| graph.blocks.iter().map(|block| &block.sequence));

  let minimizers = sequences
    .enumerate()
    .flat_map(|(i, seq)| minimizers_sketch(seq, i as u64, params))
    .sorted_by_key(|minimizer| minimizer.value)
    .collect_vec();

  let n = graphs.len();

  let mut distance = Array::zeros((n, n));
  let mut l: usize = 0;
  let mut r: usize = 0;
  while l < minimizers.len() {
    while r < minimizers.len() && minimizers[r].value == minimizers[l].value {
      r += 1;
    }

    let hits = minimizers[l..(r - 1)]
      .iter()
      .map(|m| (m.position >> 32) as usize)
      .unique()
      .sorted();

    for (i, j) in hits.tuples() {
      distance[(i, j)] += 1;
    }

    l = r;
  }

  for i in 0..n {
    for j in (i + 1)..n {
      distance[(i, j)] = 1 - distance[(i, j)] / distance[(i, i)];
      distance[(j, i)] = distance[(i, j)];
    }
    distance[(i, i)] = 0;
  }

  distance
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::pangraph::PangraphBlock;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  fn create_fake_graph(seq: impl AsRef<str>) -> Pangraph {
    let block = PangraphBlock {
      id: 0,
      sequence: seq.as_ref().to_owned(),
      gaps: BTreeMap::default(),
      mutate: vec![],
      insert: vec![],
      delete: vec![],
      positions: vec![],
    };

    let path = PangraphPath {
      name: "".to_owned(),
      nodes: vec![],
      offset: None,
      circular: false,
      position: vec![],
    };

    Pangraph {
      blocks: vec![block],
      paths: vec![path],
    }
  }

  #[rstest]
  fn test_mash_distance_general_case() {
    #[rustfmt::skip]

    let graphs = [
      create_fake_graph("ATGCATGC"),
      create_fake_graph("ATGCATGC"),
      create_fake_graph("ATGCATGC"),
      create_fake_graph("ATGCATGC"),
      create_fake_graph("ATGCATGC"),
    ];

    let actual = mash_distance(&graphs, &MinimizersParams::default());

    #[rustfmt::skip]
    let expected = array![
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
    ];

    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_mash_distance_empty() {
    let graphs = [];
    let actual = mash_distance(&graphs, &MinimizersParams::default());
    let expected = Array2::<u64>::default((0, 0));
    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_mash_distance_one() {
    let graphs = [create_fake_graph("ATGCATGC")];
    let actual = mash_distance(&graphs, &MinimizersParams::default());
    let expected = Array2::<u64>::default((0, 0));
    assert_eq!(actual, expected);
  }
}
