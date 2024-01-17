use crate::distance::mash::hash::hash;
use crate::utils::number_min_max::IsMinMaxValueOfType;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Clone, Serialize, Deserialize, Debug, SmartDefault)]
pub struct MinimizersParams {
  #[default = 15]
  pub k: u64,
  #[default = 100]
  pub w: u64,
}

/// A minimizer is a k-mer that, given a hash function that maps k-mers to integers, is the minimum k-mer within a given set of k-mers.
/// The value is the result of applying the hash function to the k-mer.
/// The position is a bit packed integer that includes reference ID, locus, and strand.
#[derive(Clone, Serialize, Deserialize, Debug, Eq, PartialEq)]
pub struct Minimizer {
  pub value: u64, // TODO: decide on a smaller the integer type here and in related functions, if memory consumption is of significance
  pub position: u64,
}

impl Minimizer {
  #[must_use]
  pub fn new(value: u64, position: u64) -> Self {
    Self { value, position }
  }

  #[must_use]
  pub fn max() -> Self {
    Self {
      value: u64::MAX,
      position: u64::MAX,
    }
  }

  pub fn is_max(&self) -> bool {
    self.value.is_max_value_of_type()
  }
}

/// Sketch a linear sequence into a vector of minimizers.
///  - `k` sets the k-mer size.
///  - `w` sets the number of contiguous k-mers that will be used in the window minimizer comparison.
///  - `id` is a unique integer that corresponds to the sequence. It will be bit-packed into the minimizer position.
pub fn minimizers_sketch(seq: impl AsRef<str>, id: u64, params: &MinimizersParams) -> Vec<Minimizer> {
  let MinimizersParams { k, w } = params;

  assert!(k < &32);
  assert!(w < &256);

  let mut fwd: u64 = 0;
  let mut rev: u64 = 0;

  let mask: u64 = (1 << (2 * k)) - 1;
  let shift: u64 = 2 * (k - 1);

  let mut min = Minimizer::max();
  let mut minimizer: Vec<Minimizer> = vec![];
  let mut window = vec![Minimizer::max(); *w as usize];

  let mut l: u64 = 0;
  let mut bi: usize = 0;
  let mut mi: usize = 0;

  for (locus, nuc) in seq.as_ref().chars().enumerate() {
    let locus = (locus as u64) + 1;
    let c = MAP[nuc as usize];

    let new = if c >= 4 {
      l = 0;
      Minimizer::max()
    } else {
      fwd = ((fwd << 2) | c) & mask;
      rev = (rev >> 2) | ((3 ^ c) << shift);
      l += 1;
      if l >= *k {
        let pos = (id << 32) | (locus << 1);
        if fwd <= rev {
          Minimizer::new(hash(fwd, mask), pos)
        } else {
          Minimizer::new(hash(rev, mask), pos | 1)
        }
      } else {
        Minimizer::max()
      }
    };

    window[bi] = new.clone();
    if (l == w + k - 1) && !min.is_max() {
      for i in (bi + 1)..(*w as usize) {
        if min.value == window[i].value && min.position != window[i].position {
          minimizer.push(window[i].clone());
        }
      }
      for i in 0..=bi {
        if min.value == window[i].value && min.position != window[i].position {
          minimizer.push(window[i].clone());
        }
      }
    }

    if new.value < min.value {
      if l >= w + k && !min.is_max() {
        minimizer.push(min);
      }
      min = new;
      mi = bi;
    } else if bi == mi {
      if (l >= w + k - 1) && !min.is_max() {
        minimizer.push(min.clone());
      }

      min = Minimizer::new(u64::MAX, min.position);
      for i in (bi + 1)..(*w as usize) {
        if window[i].value < min.value {
          mi = i;
          min = window[i].clone();
        }
      }

      for i in 0..=bi {
        if window[i].value < min.value {
          mi = i;
          min = window[i].clone();
        }
      }

      if (l >= w + k - 1) && !min.is_max() {
        for i in (bi + 1)..(*w as usize) {
          if min.value == window[i].value && min.position != window[i].position {
            minimizer.push(window[i].clone());
          }
        }
        for i in 0..=bi {
          if min.value == window[i].value && min.position != window[i].position {
            minimizer.push(window[i].clone());
          }
        }
      }
    }

    bi += 1;
    if bi >= *w as usize {
      bi = 0;
    }
  }

  if !min.is_max() {
    minimizer.push(min);
  }

  minimizer
}

#[rustfmt::skip]
const MAP: [u64; 256] = [
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_minimizers_sketch_general_case() {
    let params = MinimizersParams { w: 16, k: 8 };
    let seq = "CGATCCTTCGGGAACGTGTGACGCGAAGGTGCATGGGAGATCTCGCATTGCTGTTCTGGACGACGCGAAGAGTACTGCTACTTTCATGTCGCCTACGCCT";
    let actual = minimizers_sketch(seq, 1, &params);
    let expected = vec![
      (9685, 4294967328),
      (7669, 4294967355),
      (5583, 4294967359),
      (3600, 4294967386),
      (2383, 4294967415),
      (4791, 4294967427),
      (5338, 4294967451),
      (2190, 4294967461),
      (378, 4294967466),
    ];
    let expected = expected
      .into_iter()
      .map(|(value, position)| Minimizer { value, position })
      .collect::<Vec<_>>();
    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_minimizers_sketch_empty() {
    let actual = minimizers_sketch("", 0, &MinimizersParams::default());
    assert_eq!(actual, vec![]);
  }
}
