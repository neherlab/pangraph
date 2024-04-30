use crate::distance::mash::minimizer::{minimizers_sketch, MinimizersParams};
use crate::pangraph::pangraph::Pangraph;
use itertools::Itertools;
use ndarray::{array, Array, Array2};
use serde::{Deserialize, Serialize};

/// Compute the pairwise distance between all input graphs.
/// Distance is the set distance between minimizers.
/// Linear-time algorithm using hash collisions.
pub fn mash_distance(graphs: &[Pangraph], params: &MinimizersParams) -> Array2<f64> {
  if graphs.is_empty() {
    return array![[]];
  }

  let sequences = graphs
    .iter()
    .flat_map(|graph| graph.blocks.iter().map(|block| &block.consensus));

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

    let hits: Vec<_> = minimizers[l..r]
      .iter()
      .map(|m| (m.position >> 32) as usize)
      .unique()
      .sorted()
      .collect();

    for i in 0..hits.len() {
      for j in i..hits.len() {
        distance[(hits[i], hits[j])] += 1.0;
      }
    }

    l = r;
  }

  for i in 0..n {
    for j in (i + 1)..n {
      distance[(i, j)] = 1.0 - distance[(i, j)] / distance[(i, i)];
      distance[(j, i)] = distance[(i, j)];
    }
    distance[(i, i)] = 0.0;
  }

  distance
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::pangraph::pangraph::PangraphPath;
  use crate::pangraph::pangraph_block::PangraphBlock;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  fn create_fake_graph(seq: impl AsRef<str>) -> Pangraph {
    let block = PangraphBlock::from_consensus(seq);
    let path = PangraphPath::new(o!(""), block.id, false);
    Pangraph {
      blocks: vec![block],
      paths: vec![path],
    }
  }

  #[rstest]
  fn test_mash_distance_general_case() {
    #[rustfmt::skip]

    let params = MinimizersParams { w: 16, k: 8 };

    // from this tree:
    //    |---------1
    //    |         |--2
    // ---|            |--3
    //    |
    //    |---------4
    //              |--5
    //                 |--6

    let graphs = vec![
      "CATAGAAGCAGTCCCTGAGCACGACGCGTGTAACAATCGTTTTCAGACCTAGGACGTTAGAATATCGATCGCACGCTACGACCGACGATTAGCCGCACGAGCAAGTCGAAAACCCGAGTTAAGAGGCTGGACGTGATCCTAGACTTCGTC",
      "CATAGAAGCAGTCCCTGAGCACGAGGCGCGCAACAATCGTTTTCAGCCCTAGGACGTTAGAATATTGATCACAAGCTACGACCGACGATTAGCCGCACGAGCAAGTCGACAACCCGAGTTAAGAGGCTGGACGTGATGCTAGACTTCGTC",
      "CATAGAAGCAGTCCCTGAGCATGACGCGCGCAACGATCGTTTTCAGCCCTAGCACGTGAGAATATTGATCACAAGCTACGACCGACGATTAGCCGCACGAGCTAGTCGCCAACCCGAGTAAGGAGGCTGGACGTGATGCTAGACTACGTC",
      "ACATCAAAACTTAAAGTCGGTTACCATCTACAAATGTAGTAAGGGGGATTCTAATGAGAGAAGTGGACTGTGTAGATGGACCCGCTCACCTGCCCAGTATCTTAGTGGCGTATTCAGGATCTGGGAGGATTTGTTATTGCCTATTAGAGA",
      "ACATCAAAACTTAAAGTCGGTTCCCATCTACAAAAGTAGAAAGGGGGATTCTAATGAGAGATGTGGACTGTGTAGATGGACCCGCTAACCTGGCCAGTTTCTTAGTGGCTTAATCAGGATCTGGGAGGATTCGTTACTGCCTATTAGAGA",
      "ACATCAGAACTTAAAGTCGGTTCCTATCTCCAAAAGTATAAAGTGGGATTCTAATGAGAGATGTGGACTGTGTCGATAAACCCGCTAACCTGGCCTGTTTCTTGTTGGCTTAATCAGGATCTGAGAGGATTCGTTACTGCCTAGTAGTGA",
    ].into_iter().map(create_fake_graph).collect_vec();

    let actual = mash_distance(&graphs, &params);

    let expected = array![
      [0.0, 1. - 6. / 9., 0.75, 1.0, 1.0, 1.0],
      [1. - 6. / 9., 0.0, 0.5, 1.0, 1.0, 1.0],
      [0.75, 0.5, 0.0, 1.0, 1.0, 1.0],
      [1.0, 1.0, 1.0, 0.0, 0.625, 0.875],
      [1.0, 1.0, 1.0, 0.625, 0.0, 5. / 7.],
      [1.0, 1.0, 1.0, 0.875, 5. / 7., 0.0],
    ];

    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_mash_distance_empty() {
    let graphs = [];
    let actual = mash_distance(&graphs, &MinimizersParams::default());
    let expected: Array2<f64> = array![[]];
    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_mash_distance_one() {
    let graphs = [create_fake_graph("ATGCATGC")];
    let actual = mash_distance(&graphs, &MinimizersParams::default());
    let expected = array![[0.0]];
    assert_eq!(actual, expected);
  }
}
