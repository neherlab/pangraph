/// A transliteration of Jenkin's invertible hash function for 64 bit integers.
/// Bijectively maps any k-mer to an integer.
pub const fn hash(x: u64, mask: u64) -> u64 {
  let mut x = (!x).wrapping_add(x << 21) & mask;
  x = x ^ (x >> 24);
  x = (x + (x << 3) + (x << 8)) & mask;
  x = x ^ (x >> 14);
  x = (x + (x << 2) + (x << 4)) & mask;
  x = x ^ (x >> 28);
  x = (x + (x << 31)) & mask;
  x
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  #[case::zero((0, 0), 0)]
  #[case::zero_mask((123, 0), 0)]
  #[case::zero_value((0, 456), 136)]
  #[case::random((123, 456), 384)]
  fn test_hash(#[case] input: (u64, u64), #[case] expected: u64) {
    assert_eq!(hash(input.0, input.1), expected);
  }
}
