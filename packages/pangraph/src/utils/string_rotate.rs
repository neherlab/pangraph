pub trait StringRotateLeft {
  fn rotate_left(&mut self, mid: usize);
}

impl StringRotateLeft for String {
  fn rotate_left(&mut self, mid: usize) {
    assert!(mid <= self.len());
    debug_assert!(self.is_ascii());

    // SAFETY: String must be valid UTF-8 when string is used.
    // In our case the char set is ASCII, so it does not change after rotation.
    #[allow(unsafe_code)]
    let slice = unsafe { self.as_bytes_mut() };
    slice.rotate_left(mid);
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{o, pretty_assert_eq};
  use rstest::rstest;

  #[rstest]
  #[case("hello world", 5, " worldhello")]
  #[case("hello world", 3, "lo worldhel")]
  #[case("hello world", 0, "hello world")]
  #[case("hello world", 11, "hello world")] // Full Rotation
  #[case("", 0, "")]
  fn test_string_rotate_left_basic(#[case] mut input: String, #[case] mid: usize, #[case] expected: &str) {
    input.rotate_left(mid);
    pretty_assert_eq!(input, expected);
  }

  #[test]
  fn test_string_rotate_left_idempotence() {
    let mut input = o!("hello world");
    let original = input.clone();
    input.rotate_left(input.len());
    pretty_assert_eq!(input, original);
  }

  #[test]
  fn test_string_rotate_left_associativity() {
    let mut input = o!("hello world");
    let original = input.clone();
    input.rotate_left(3);
    input.rotate_left(5);
    let mut expected = original;
    expected.rotate_left(8); // Combine rotations (3 + 5 = 8)
    assert_eq!(input, expected);
  }
}
