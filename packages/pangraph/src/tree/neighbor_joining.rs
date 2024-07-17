#![allow(non_snake_case)]

use crate::commands::build::build_args::PangraphBuildArgs;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GetName, GraphNode, GraphNodeKey};
use crate::pangraph::pangraph::Pangraph;
use crate::tree::distances::{calculate_distances, create_Q_matrix, dist, pair};
use derive_more::Display;
use eyre::Report;
use itertools::Itertools;
use ndarray::{array, s, Array2, AssignElem, Axis};
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};
use std::ops::Deref;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Clade {
  pub pangraph: Option<Pangraph>,
}

impl GraphNode for Clade {}

#[derive(Clone, Debug, Default, Display, Serialize, Deserialize)]
pub struct Empty;

impl GraphEdge for Empty {}

pub type Tree = Graph<Clade, Empty>;

/// Generate guide tree using neighbor joining method.
pub fn build_tree_using_neighbor_joining(graphs: Vec<Pangraph>, args: &PangraphBuildArgs) -> Result<Tree, Report> {
  let mut distances = calculate_distances(&graphs, args);

  let mut tree = Tree::new();

  let mut node_keys = graphs
    .into_iter()
    .map(|pangraph| {
      tree.add_node(Clade {
        pangraph: Some(pangraph),
      })
    })
    .collect_vec();

  join_all_in_place(&mut tree, &mut node_keys, &mut distances)?;

  // Balance guide tree (to increase available parallelism during parallel traversal?)
  // let tree = balance(&tree);

  Ok(tree)
}

fn join_all_in_place<N: GraphNode + Default, E: GraphEdge + Default>(
  tree: &mut Graph<N, E>,
  node_keys: &mut Vec<GraphNodeKey>,
  distances: &mut Array2<f64>,
) -> Result<(), Report> {
  while distances.len() > 4 {
    join_in_place(tree, node_keys, distances)?;
  }

  let root = tree.add_node(N::default());
  tree.add_edge(root, node_keys[0], E::default())?;
  tree.add_edge(root, node_keys[1], E::default())?;

  tree.build()?;
  Ok(())
}

fn join_in_place<N: GraphNode + Default, E: GraphEdge + Default>(
  tree: &mut Graph<N, E>,
  node_keys: &mut Vec<GraphNodeKey>,
  D: &mut Array2<f64>,
) -> Result<(), Report> {
  let q = create_Q_matrix(D).unwrap();
  let (i, j) = pair(&q)?;

  let node_key = tree.add_node(N::default());
  tree.add_edge(node_key, node_keys[i], E::default())?;
  tree.add_edge(node_key, node_keys[j], E::default())?;
  node_keys[i] = node_key;
  node_keys.remove(j);

  let dn = dist(D, i, j);
  D.slice_mut(s![i, ..]).assign(&dn);
  D.slice_mut(s![.., i]).assign(&dn);
  D[(i, i)] = 0.0;

  D.remove_index(Axis(0), j);
  D.remove_index(Axis(1), j);

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::io::nwk::tests::{TestEdge, TestNode};
  use crate::io::nwk::{nwk_write_str, WriteNwkOptions};
  use ndarray::array;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_tree_neighbor_joining_distances() -> Result<(), Report> {
    let mut tree = Graph::<TestNode, TestEdge>::new();

    let mut nodes = vec!["A", "B", "C", "D", "E"]
      .into_iter()
      .map(|name| tree.add_node(TestNode(name.to_owned())))
      .collect_vec();

    #[rustfmt::skip]
    let mut D = array![
      [0.0,  5.0,  9.0,  9.0, 8.0],
      [5.0,  0.0, 10.0, 10.0, 9.0],
      [9.0, 10.0,  0.0,  8.0, 7.0],
      [9.0, 10.0,  8.0,  0.0, 3.0],
      [8.0,  9.0,  7.0,  3.0, 0.0],
    ];

    join_in_place(&mut tree, &mut nodes, &mut D)?;

    #[rustfmt::skip]
    let D_expected = array![
      [0.0, 7.0,  7.0, 6.0],
      [7.0, 0.0,  8.0, 7.0],
      [7.0, 8.0,  0.0, 3.0],
      [6.0, 7.0,  3.0, 0.0],
    ];

    assert_eq!(&D_expected, &D);

    join_in_place(&mut tree, &mut nodes, &mut D)?;

    #[rustfmt::skip]
    let D_expected = array![
      [0.0,  4.0, 3.0],
      [4.0,  0.0, 3.0],
      [3.0,  3.0, 0.0],
    ];

    assert_eq!(&D, &D_expected);

    Ok(())
  }

  #[test]
  fn test_tree_neighbor_joining_tree_simple() -> Result<(), Report> {
    let mut tree = Graph::<TestNode, TestEdge>::new();

    let mut nodes = vec!["A", "B", "C", "D", "E"]
      .into_iter()
      .map(|name| tree.add_node(TestNode(name.to_owned())))
      .collect_vec();

    #[rustfmt::skip]
    let mut distances = array![
      [0.0,  5.0,  9.0,  9.0, 8.0],
      [5.0,  0.0, 10.0, 10.0, 9.0],
      [9.0, 10.0,  0.0,  8.0, 7.0],
      [9.0, 10.0,  8.0,  0.0, 3.0],
      [8.0,  9.0,  7.0,  3.0, 0.0],
    ];

    join_all_in_place(&mut tree, &mut nodes, &mut distances)?;

    let actual = nwk_write_str(&tree, &WriteNwkOptions::default())?;

    assert_eq!("((((A,B),C),D),E);", actual);

    Ok(())
  }

  #[test]
  fn test_tree_neighbor_joining_tree_general() -> Result<(), Report> {
    let mut tree = Graph::<TestNode, TestEdge>::new();

    let mut nodes = vec!["A", "B", "C", "D", "E", "F", "G", "H"]
      .into_iter()
      .map(|name| tree.add_node(TestNode(name.to_owned())))
      .collect_vec();

    #[rustfmt::skip]
    let mut distances = array![
      //  A     B     C     D     E     F     G     H
      [ 0.0, 46.0, 37.0, 46.0, 46.0, 14.0, 37.0,  1.0],  // A
      [46.0,  0.0, 46.0,  7.0,  1.0, 46.0, 46.0, 46.0],  // B
      [37.0, 46.0,  0.0, 46.0, 46.0, 37.0,  1.0, 37.0],  // C
      [46.0,  7.0, 46.0,  0.0,  7.0, 46.0, 46.0, 46.0],  // D
      [46.0,  1.0, 46.0,  7.0,  0.0, 46.0, 46.0, 46.0],  // E
      [14.0, 46.0, 37.0, 46.0, 46.0,  0.0, 37.0, 14.0],  // F
      [37.0, 46.0,  1.0, 46.0, 46.0, 37.0,  0.0, 37.0],  // G
      [ 1.0, 46.0, 37.0, 46.0, 46.0, 14.0, 37.0,  0.0],  // H
    ];

    join_all_in_place(&mut tree, &mut nodes, &mut distances)?;

    let actual = nwk_write_str(&tree, &WriteNwkOptions::default())?;

    assert_eq!("(((A,H),(((B,E),D),(C,G))),F);", actual);

    // expected tree:
    //               __ A
    //      ________|
    //     |        |__ H
    //     |
    //     |          __ B
    //     |       __|
    //   __|      |  |__ E
    //  |  |   ___|
    //  |  |  |   |__ D
    //  |  |__|
    //  |     |    __ C
    //  |     |___|
    //  |         |__ G
    //  |
    //  |____________ F
    //
    // newick: (((A,H),(((B,E),D),(C,G))),F)

    Ok(())
  }
}
