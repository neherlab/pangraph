#![allow(non_snake_case)]

use crate::commands::build::build_args::{DistanceBackend, PangraphBuildArgs};
use crate::distance::mash::mash_distance::mash_distance;
use crate::distance::mash::minimizer::MinimizersParams;
use crate::pangraph::pangraph::Pangraph;
use crate::utils::ndarray::broadcast;
use eyre::Report;
use ndarray::{s, Array1, Array2, Axis};
use ndarray_stats::QuantileExt;

// Calculate pairwise distances between future guide tree nodes
pub fn calculate_distances(graphs: &[Pangraph], args: &PangraphBuildArgs) -> Array2<f64> {
  let distances = match args.distance_backend {
    // TODO: this function only needs sequences, and not graphs
    DistanceBackend::Native => mash_distance(graphs, &MinimizersParams::default()),
    DistanceBackend::Mash => {
      // FIXME: what's the difference between Native and Mash?
      unimplemented!("DistanceBackend::Mash");
    }
  };

  assert_eq!(distances.len_of(Axis(0)), distances.len_of(Axis(1)));
  distances
}

pub fn create_Q_matrix(D: &Array2<f64>) -> Result<Array2<f64>, Report> {
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

pub fn pair(Q: &Array2<f64>) -> Result<(usize, usize), Report> {
  let iota = Q.argmin()?;
  Ok(if iota.0 > iota.1 {
    (iota.1, iota.0)
  } else {
    (iota.0, iota.1)
  })
}

pub fn dist(D: &Array2<f64>, i: usize, j: usize) -> Array1<f64> {
  let n = D.len_of(Axis(0));
  assert!(n > 2);

  let dn: Array1<f64> = 0.5 * (&D.slice(s![i, ..]) + &D.slice(s![j, ..]) - D[(i, j)]);

  dn
}

#[cfg(test)]
mod tests {
  use super::*;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  const INF: f64 = f64::INFINITY;

  #[test]
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

  #[test]
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
}
