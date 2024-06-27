#![allow(non_snake_case)]

use crate::make_error;
use crate::pangraph::pangraph::Pangraph;
use crate::tree::clade::Clade;
use crate::utils::lock::Lock;
use crate::utils::ndarray::broadcast;
use eyre::Report;
use itertools::Itertools;
use ndarray::{array, s, Array1, Array2, AssignElem, Axis};
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};

/// Generate a tree from a matrix of pairwise distances `distance` using neighbor joining method.
/// The names of leafs are given by an array of strings `names`.
pub fn build_tree_using_neighbor_joining(
  distance: &Array2<f64>,
  names: &[String],
  graphs: &[Pangraph],
) -> Result<Lock<Clade>, Report> {
  assert_eq!(names.len(), graphs.len());
  assert_eq!(names.len(), distance.len_of(Axis(0)));
  assert_eq!(distance.len_of(Axis(0)), distance.len_of(Axis(1)));

  if names.len() < 2 {
    return make_error!("Expected at least 2 samples, but found {}", names.len());
  }

  let mut nodes = names
    .iter()
    .map(|name| {
      // // TODO: try to restructure the code such that there is no need to store graphs and names in
      // // separate arrays and do fallible lookups. The relation is 1:1 and there should be no need for that.
      // // Also external names are unreliable and are not guaranteed to uniquely identify the sequences.
      // // TODO: Looking up names in the paths make testing setup of this function harder.

      // let graph = graphs
      //   .iter()
      //   .find(|graph| &graph.paths[&PathId(0)].name == name)
      //   .ok_or_else(|| make_internal_report!("Graph not found for clade '{name}'"))?;
      let graph = &graphs[0]; // FIXME

      Ok(Lock::new(Clade::new(name, Some(graph.clone()))))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let mut distance = distance.clone(); // TODO: should we avoid copying here?

  while nodes.len() > 2 {
    join_in_place(&mut distance, &mut nodes)?;
  }

  Ok(Lock::new(Clade::from_children(&nodes[0], &nodes[1])))
}

fn create_Q_matrix(D: &Array2<f64>) -> Result<Array2<f64>, Report> {
  assert_eq!(D.len_of(Axis(0)), D.len_of(Axis(1)));

  let n = D.len_of(Axis(0));
  assert!(n > 2);

  let sum_0 = D.sum_axis(Axis(0));
  let sum_0 = broadcast(&sum_0, (n, n))?;

  let sum_1 = D.sum_axis(Axis(1));
  let sum_1 = broadcast(&sum_1, (n, n))?;

  let mut Q = ((n as f64) - 2.0) * D - sum_0 - sum_1.t();
  Q.diag_mut().fill(f64::INFINITY);
  Ok(Q)
}

fn create_Q_matrix_native(D: &Array2<f64>) -> Array2<f64> {
  let n = D.len_of(Axis(0));
  assert!(n > 2);

  let mut Q = Array2::zeros((n, n));
  for i in 0..n {
    for j in 0..n {
      let mut sum_0 = 0.0;
      let mut sum_1 = 0.0;
      for k in 0..n {
        sum_0 += D[(i, k)];
        sum_1 += D[(j, k)];
      }
      Q[(i, j)] = ((n as f64) - 2.0) * D[(i, j)] - sum_0 - sum_1;
    }
  }

  Q.diag_mut().fill(f64::INFINITY);
  Q
}

fn pair(Q: &Array2<f64>) -> Result<(usize, usize), Report> {
  let iota = Q.argmin()?;
  Ok(if iota.0 > iota.1 {
    (iota.1, iota.0)
  } else {
    (iota.0, iota.1)
  })
}

fn dist(D: &Array2<f64>, i: usize, j: usize) -> Array1<f64> {
  let n = D.len_of(Axis(0));
  assert!(n > 2);

  let dn: Array1<f64> = 0.5 * (&D.slice(s![i, ..]) + &D.slice(s![j, ..]) - D[(i, j)]);

  dn
}

fn join_in_place(D: &mut Array2<f64>, nodes: &mut Vec<Lock<Clade>>) -> Result<(), Report> {
  let q = create_Q_matrix(D).unwrap();
  let (i, j) = pair(&q)?;

  let node = Lock::new(Clade::from_children(&nodes[i], &nodes[j]));
  nodes[i].write().parent = Some(node.clone()); // TODO: is this assignment redundant? (node[i] is overwritten below)
  nodes[j].write().parent = Some(node.clone());
  nodes[i] = node;
  nodes.remove(j);

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
  use crate::o;
  use itertools::Itertools;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const INF: f64 = f64::INFINITY;

  #[rstest]
  fn test_create_Q_matrix() {
    // example from wikipedia: https://en.wikipedia.org/wiki/Neighbor_joining#Example
    #[rustfmt::skip]
    let distance = array![
      [0.0,  5.0,  9.0,  9.0, 8.0],
      [5.0,  0.0, 10.0, 10.0, 9.0],
      [9.0, 10.0,  0.0,  8.0, 7.0],
      [9.0, 10.0,  8.0,  0.0, 3.0],
      [8.0,  9.0,  7.0,  3.0, 0.0],
    ];

    let actual = create_Q_matrix(&distance).unwrap();

    #[rustfmt::skip]
    let expected = array![
      [  INF, -50.0, -38.0 , -34.0, -34.0],
      [-50.0,   INF, -38.0,  -34.0, -34.0],
      [-38.0, -38.0,   INF , -40.0, -40.0],
      [-34.0, -34.0,  -40.0,   INF, -48.0],
      [-34.0, -34.0,  -40.0, -48.0,   INF],
    ];

    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_create_Q_matrix_native() {
    #[rustfmt::skip]
    let distance = array![
      [0.0,  5.0,  9.0,  9.0, 8.0],
      [5.0,  0.0, 10.0, 10.0, 9.0],
      [9.0, 10.0,  0.0,  8.0, 7.0],
      [9.0, 10.0,  8.0,  0.0, 3.0],
      [8.0,  9.0,  7.0,  3.0, 0.0],
    ];

    let actual = create_Q_matrix_native(&distance);

    #[rustfmt::skip]
    let expected = array![
      [  INF, -50.0, -38.0 , -34.0, -34.0],
      [-50.0,   INF, -38.0,  -34.0, -34.0],
      [-38.0, -38.0,   INF , -40.0, -40.0],
      [-34.0, -34.0,  -40.0,   INF, -48.0],
      [-34.0, -34.0,  -40.0, -48.0,   INF],
    ];

    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_dist() {
    #[rustfmt::skip]
    let distance = array![
      [0.0,  5.0,  9.0,  9.0, 8.0],
      [5.0,  0.0, 10.0, 10.0, 9.0],
      [9.0, 10.0,  0.0,  8.0, 7.0],
      [9.0, 10.0,  8.0,  0.0, 3.0],
      [8.0,  9.0,  7.0,  3.0, 0.0],
    ];

    let actual = dist(&distance, 0, 1);
    let expected = array![0., 0., 7., 7., 6.];
    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_join() {
    let mut nodes = vec!["A", "B", "C", "D", "E"]
      .into_iter()
      .map(|name| Lock::new(Clade::new(name, None)))
      .collect_vec();

    #[rustfmt::skip]
    let mut D = array![
      [0.0,  5.0,  9.0,  9.0, 8.0],
      [5.0,  0.0, 10.0, 10.0, 9.0],
      [9.0, 10.0,  0.0,  8.0, 7.0],
      [9.0, 10.0,  8.0,  0.0, 3.0],
      [8.0,  9.0,  7.0,  3.0, 0.0],
    ];

    join_in_place(&mut D, &mut nodes).unwrap();

    #[rustfmt::skip]
    let D_expected = array![
      [0.0, 7.0,  7.0, 6.0],
      [7.0, 0.0,  8.0, 7.0],
      [7.0, 8.0,  0.0, 3.0],
      [6.0, 7.0,  3.0, 0.0],
    ];

    assert_eq!(&D, &D_expected);

    let nodes_actual = nodes.iter().map(|node| node.read_arc().name.clone()).collect_vec();
    let nodes_expected = vec![None, Some(o!("C")), Some(o!("D")), Some(o!("E"))];
    assert_eq!(nodes_actual, nodes_expected);

    join_in_place(&mut D, &mut nodes).unwrap();

    #[rustfmt::skip]
    let D_expected = array![
      [0.0,  4.0, 3.0],
      [4.0,  0.0, 3.0],
      [3.0,  3.0, 0.0],
    ];

    assert_eq!(&D, &D_expected);

    let nodes_actual = nodes.iter().map(|node| node.read_arc().name.clone()).collect_vec();
    let nodes_expected = vec![None, Some(o!("D")), Some(o!("E"))];
    assert_eq!(nodes_actual, nodes_expected);
  }

  // FIXME
  //
  // #[rstest]
  // fn test_build_tree_using_neighbor_joining() {
  //   let names = vec_of_owned!["A", "B", "C", "D", "E"];
  //   let graphs = names
  //     .iter()
  //     .map(|name| Pangraph {
  //       paths: vec![],
  //       blocks: vec![],
  //     })
  //     .collect_vec();
  //
  //   #[rustfmt::skip]
  //   let D = array![
  //     [0.0,  5.0,  9.0,  9.0, 8.0],
  //     [5.0,  0.0, 10.0, 10.0, 9.0],
  //     [9.0, 10.0,  0.0,  8.0, 7.0],
  //     [9.0, 10.0,  8.0,  0.0, 3.0],
  //     [8.0,  9.0,  7.0,  3.0, 0.0],
  //   ];
  //
  //   let tree = build_tree_using_neighbor_joining(&D, &names, &graphs).unwrap();
  //
  //   let nwk = tree.read().to_newick();
  //   assert_eq!(nwk, "((((A,B),C),D),E);");
  //
  //   let nodes: Vec<String> = postorder(&tree, |clade| clade.name.clone().unwrap_or_default());
  //
  //   let nodes_expected: Vec<String> = vec_of_owned!["A", "B", "", "C", "", "D", "", "E", ""];
  //   assert_eq!(nodes, nodes_expected);
  // }

  // FIXME
  //
  // #[rstest]
  // fn test_build_tree_using_neighbor_joining_2() {
  //   let names = vec_of_owned!["A", "B", "C", "D", "E", "F", "G", "H"];
  //   let graphs = names
  //     .iter()
  //     .map(|name| Pangraph {
  //       paths: vec![],
  //       blocks: vec![],
  //     })
  //     .collect_vec();
  //
  //   // expected tree:
  //   //               __ A
  //   //      ________|
  //   //     |        |__ H
  //   //     |
  //   //     |          __ B
  //   //     |       __|
  //   //   __|      |  |__ E
  //   //  |  |   ___|
  //   //  |  |  |   |__ D
  //   //  |  |__|
  //   //  |     |    __ C
  //   //  |     |___|
  //   //  |         |__ G
  //   //  |
  //   //  |____________ F
  //   //
  //   // newick: (((A,H),(((B,E),D),(C,G))),F)
  //
  //   #[rustfmt::skip]
  //   let D = array![
  //     //  A     B     C     D     E     F     G     H
  //     [ 0.0, 46.0, 37.0, 46.0, 46.0, 14.0, 37.0,  1.0],  // A
  //     [46.0,  0.0, 46.0,  7.0,  1.0, 46.0, 46.0, 46.0],  // B
  //     [37.0, 46.0,  0.0, 46.0, 46.0, 37.0,  1.0, 37.0],  // C
  //     [46.0,  7.0, 46.0,  0.0,  7.0, 46.0, 46.0, 46.0],  // D
  //     [46.0,  1.0, 46.0,  7.0,  0.0, 46.0, 46.0, 46.0],  // E
  //     [14.0, 46.0, 37.0, 46.0, 46.0,  0.0, 37.0, 14.0],  // F
  //     [37.0, 46.0,  1.0, 46.0, 46.0, 37.0,  0.0, 37.0],  // G
  //     [ 1.0, 46.0, 37.0, 46.0, 46.0, 14.0, 37.0,  0.0],  // H
  //   ];
  //
  //   let tree = build_tree_using_neighbor_joining(&D, &names, &graphs).unwrap();
  //
  //   let nwk = tree.read().to_newick();
  //   assert_eq!(nwk, "(((A,H),(((B,E),D),(C,G))),F);");
  //
  //   let nodes: Vec<String> = postorder(&tree, |clade| clade.name.clone().unwrap_or_default());
  //
  //   let nodes_expected: Vec<String> = vec_of_owned!["A", "H", "", "B", "E", "", "D", "", "C", "G", "", "", "", "F", ""];
  //   assert_eq!(nodes, nodes_expected);
  // }
}
