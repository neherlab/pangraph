use itertools::Itertools;
use std::collections::BTreeMap;

pub type Fun<K, V> = fn((K, V), (K, V)) -> V;

#[derive(Copy, Clone)]
pub enum ConflictResolution<K, V> {
  Left,
  Right,
  Custom(Fun<K, V>),
}

#[allow(clippy::needless_pass_by_value)]
#[must_use]
pub fn map_merge<K, V>(
  left: &BTreeMap<K, V>,
  right: &BTreeMap<K, V>,
  resolution: ConflictResolution<K, V>,
) -> BTreeMap<K, V>
where
  K: Ord + Clone,
  V: Clone,
{
  left
    .iter()
    .merge_join_by(right.iter(), |(kl, _), (kr, _)| kl.cmp(kr))
    .map(|either| match either {
      itertools::EitherOrBoth::Left((k, v)) => (k.clone(), v.clone()),
      itertools::EitherOrBoth::Right((k, v)) => (k.clone(), v.clone()),
      itertools::EitherOrBoth::Both((kl, vl), (kr, vr)) => match resolution {
        ConflictResolution::Left => (kl.clone(), vl.clone()),
        ConflictResolution::Right => (kr.clone(), vr.clone()),
        ConflictResolution::Custom(f) => (kl.clone(), f((kl.clone(), vl.clone()), (kr.clone(), vr.clone()))),
      },
    })
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use maplit::btreemap;

  #[test]
  fn test_map_merge_both_empty() {
    let merged = map_merge(&btreemap! {}, &btreemap! {}, ConflictResolution::<i32, i32>::Left);
    assert_eq!(merged, btreemap! {});
  }

  #[test]
  fn test_map_merge_left_map_empty_left() {
    let merged = map_merge(&btreemap! {1 => 10, 2 => 20}, &btreemap! {}, ConflictResolution::Left);
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20});
  }

  #[test]
  fn test_map_merge_right_map_empty_left() {
    let merged = map_merge(&btreemap! {}, &btreemap! {1 => 10, 2 => 20}, ConflictResolution::Left);
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20});
  }

  #[test]
  fn test_map_merge_left_map_empty_right() {
    let merged = map_merge(&btreemap! {}, &btreemap! {1 => 10, 2 => 20}, ConflictResolution::Right);
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20});
  }

  #[test]
  fn test_map_merge_right_map_empty_right() {
    let merged = map_merge(&btreemap! {1 => 10, 2 => 20}, &btreemap! {}, ConflictResolution::Right);
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20});
  }

  #[test]
  fn test_map_merge_non_overlapping_maps_left() {
    let merged = map_merge(
      &btreemap! {1 => 10, 3 => 30},
      &btreemap! {2 => 20, 4 => 40},
      ConflictResolution::Left,
    );
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20, 3 => 30, 4 => 40});
  }

  #[test]
  fn test_map_merge_non_overlapping_maps_right() {
    let merged = map_merge(
      &btreemap! {1 => 10, 3 => 30},
      &btreemap! {2 => 20, 4 => 40},
      ConflictResolution::Right,
    );
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20, 3 => 30, 4 => 40});
  }

  #[test]
  fn test_map_merge_one_key_overlap_left() {
    let merged = map_merge(
      &btreemap! {1 => 100, 2 => 200},
      &btreemap! {2 => 300, 3 => 400},
      ConflictResolution::Left,
    );
    assert_eq!(merged, btreemap! {1 => 100, 2 => 200, 3 => 400});
  }

  #[test]
  fn test_map_merge_one_key_overlap_right() {
    let merged = map_merge(
      &btreemap! {1 => 100, 2 => 200},
      &btreemap! {2 => 300, 3 => 400},
      ConflictResolution::Right,
    );
    assert_eq!(merged, btreemap! {1 => 100, 2 => 300, 3 => 400});
  }

  #[test]
  fn test_map_merge_multiple_keys_left() {
    let merged = map_merge(
      &btreemap! {1 => 10, 4 => 40, 6 => 60},
      &btreemap! {2 => 20, 4 => 45, 5 => 50, 6 => 65},
      ConflictResolution::Left,
    );
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20, 4 => 40, 5 => 50, 6 => 60});
  }

  #[test]
  fn test_map_merge_multiple_keys_right() {
    let merged = map_merge(
      &btreemap! {1 => 10, 4 => 40, 6 => 60},
      &btreemap! {2 => 20, 4 => 45, 5 => 50, 6 => 65},
      ConflictResolution::Right,
    );
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20, 4 => 45, 5 => 50, 6 => 65});
  }

  #[test]
  fn test_map_merge_multiple_keys_custom() {
    let merged = map_merge(
      &btreemap! {1 => 10, 4 => 40, 6 => 60},
      &btreemap! {2 => 20, 4 => 45, 5 => 50, 6 => 65},
      ConflictResolution::Custom(|(_, l), (_, r)| l + r),
    );
    assert_eq!(merged, btreemap! {1 => 10, 2 => 20, 4 => 85, 5 => 50, 6 => 125});
  }
}
