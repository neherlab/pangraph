use crate::representation::seq_char::AsciiChar;
use schemars::{JsonSchema, Schema, SchemaGenerator};
use std::borrow::Cow;

/// Represents genetic sequence (ASCII characters only)
#[must_use]
#[derive(Clone, PartialOrd, Ord, PartialEq, Eq, Hash, Default)]
pub struct Seq {
  data: Vec<AsciiChar>,
}

impl Seq {
  pub fn new() -> Self {
    Self { data: Vec::new() }
  }

  pub fn with_capacity(capacity: usize) -> Self {
    Self {
      data: Vec::with_capacity(capacity),
    }
  }

  pub fn from_elem<T: Into<u8>>(elem: T, n: usize) -> Self {
    Self {
      data: vec![AsciiChar::from(elem.into()); n],
    }
  }

  pub fn from_str(s: &str) -> Self {
    debug_assert!(s.is_ascii());
    Self {
      data: s.as_bytes().iter().copied().map(AsciiChar::from).collect(),
    }
  }

  pub fn from_vec(vec: Vec<u8>) -> Self {
    Self {
      data: vec.into_iter().map(AsciiChar::from).collect(),
    }
  }

  pub fn from_slice(slice: &[u8]) -> Self {
    Self {
      data: slice.iter().copied().map(AsciiChar::from).collect(),
    }
  }

  pub fn concat(sequences: &[&Self]) -> Self {
    let total_len = sequences.iter().map(|s| s.len()).sum();
    let mut data = Vec::with_capacity(total_len);
    for seq in sequences {
      data.extend_from_slice(&seq.data);
    }
    Seq { data }
  }

  pub fn extend_seq(&mut self, other: &Seq) {
    self.data.extend_from_slice(&other.data);
  }

  pub fn len(&self) -> usize {
    self.data.len()
  }

  pub fn is_empty(&self) -> bool {
    self.data.is_empty()
  }

  pub fn clear(&mut self) {
    self.data.clear();
  }

  pub fn truncate(&mut self, len: usize) {
    self.data.truncate(len);
  }

  pub fn push(&mut self, byte: AsciiChar) {
    self.data.push(byte);
  }

  pub fn pop(&mut self) -> Option<AsciiChar> {
    self.data.pop().map(Into::into)
  }

  pub fn capacity(&self) -> usize {
    self.data.capacity()
  }

  pub fn reserve(&mut self, additional: usize) {
    self.data.reserve(additional);
  }

  pub fn reserve_exact(&mut self, additional: usize) {
    self.data.reserve_exact(additional);
  }

  #[allow(unsafe_code)]
  pub fn as_str(&self) -> &str {
    // SAFETY: `self.data.as_ptr()` is guaranteed to be valid for reads and properly aligned
    // because `data` is a `Vec<AsciiChar>`. The `AsciiChar` type ensures that each element is a valid
    // single-byte ASCII character, making the conversion to a `u8` pointer valid.
    let byte_slice = unsafe { std::slice::from_raw_parts(self.data.as_ptr().cast::<u8>(), self.data.len()) };

    // SAFETY: `from_utf8_unchecked` is safe here because `byte_slice` is guaranteed to contain only
    // valid UTF-8 data. This is ensured by the invariant that `AsciiChar` can only hold valid ASCII characters,
    // which are a subset of UTF-8.
    unsafe { std::str::from_utf8_unchecked(byte_slice) }
  }

  pub fn as_slice(&self) -> &[AsciiChar] {
    &self.data
  }

  pub fn as_mut_slice(&mut self) -> &mut [AsciiChar] {
    &mut self.data
  }

  pub fn insert(&mut self, index: usize, byte: u8) {
    self.data.insert(index, AsciiChar::from(byte));
  }

  pub fn insert_seq(&mut self, index: usize, other: &Seq) {
    let len = self.len();
    assert!(
      index <= len,
      "Attempted to insert outside of array boundaries: array size is {len}, index is {index}"
    );

    self.data.reserve(other.data.len());
    let mut tail = self.data.split_off(index);
    self.data.extend_from_slice(&other.data);
    self.data.append(&mut tail);
  }

  pub fn remove(&mut self, index: usize) -> u8 {
    self.data.remove(index).into()
  }

  pub fn append(&mut self, other: &mut Vec<u8>) {
    self.data.extend(other.drain(..).map(AsciiChar::from));
  }

  pub fn push_str(&mut self, s: &str) {
    self.data.extend(s.as_bytes().iter().copied().map(AsciiChar::from));
  }

  pub fn contains_str(&self, s: &str) -> bool {
    self.as_str().contains(s)
  }

  pub fn starts_with(&self, prefix: &str) -> bool {
    self.as_str().starts_with(prefix)
  }

  pub fn ends_with(&self, suffix: &str) -> bool {
    self.as_str().ends_with(suffix)
  }

  pub fn find(&self, substring: &str) -> Option<usize> {
    self.as_str().find(substring)
  }

  pub fn rfind(&self, substring: &str) -> Option<usize> {
    self.as_str().rfind(substring)
  }

  pub fn split_off(&mut self, at: usize) -> Seq {
    Seq {
      data: self.data.split_off(at),
    }
  }

  pub fn retain<F>(&mut self, f: F)
  where
    F: FnMut(&AsciiChar) -> bool,
  {
    self.data.retain(f);
  }

  pub fn rotate_left(&mut self, mid: usize) {
    self.data.rotate_left(mid);
  }

  pub fn rotate_right(&mut self, mid: usize) {
    self.data.rotate_right(mid);
  }

  pub fn n_ambiguous_bases(&self) -> usize {
    self.data.iter().filter(|&c| c.is_ambiguous()).count()
  }

  pub fn ambiguous_positions(&self) -> Vec<usize> {
    self
      .data
      .iter()
      .enumerate()
      .filter_map(|(i, c)| c.is_ambiguous().then_some(i))
      .collect()
  }
}

impl PartialEq<str> for Seq {
  fn eq(&self, other: &str) -> bool {
    if self.data.len() != other.len() {
      return false;
    }
    self.data.iter().zip(other.as_bytes()).all(|(&AsciiChar(c), &b)| c == b)
  }
}

impl PartialEq<String> for Seq {
  fn eq(&self, other: &String) -> bool {
    self == other.as_str()
  }
}

impl PartialEq<&str> for Seq {
  fn eq(&self, other: &&str) -> bool {
    self == *other
  }
}

impl PartialEq<Seq> for str {
  fn eq(&self, other: &Seq) -> bool {
    other == self
  }
}

impl PartialEq<Seq> for String {
  fn eq(&self, other: &Seq) -> bool {
    other == self.as_str()
  }
}

impl PartialEq<&String> for Seq {
  fn eq(&self, other: &&String) -> bool {
    self == other.as_str()
  }
}

impl PartialEq<&Seq> for Seq {
  fn eq(&self, other: &&Seq) -> bool {
    self == *other
  }
}

impl PartialEq<Seq> for &str {
  fn eq(&self, other: &Seq) -> bool {
    *self == other
  }
}

impl PartialEq<Seq> for &String {
  fn eq(&self, other: &Seq) -> bool {
    self.as_str() == other
  }
}

impl PartialEq<&Seq> for str {
  fn eq(&self, other: &&Seq) -> bool {
    *other == self
  }
}

impl PartialEq<&Seq> for String {
  fn eq(&self, other: &&Seq) -> bool {
    self.as_str() == *other
  }
}

impl PartialEq<[char]> for Seq {
  fn eq(&self, other: &[char]) -> bool {
    self.data.len() == other.len() && self.data.iter().zip(other).all(|(&AsciiChar(c), &ch)| c as char == ch)
  }
}

impl PartialEq<&[char]> for Seq {
  fn eq(&self, other: &&[char]) -> bool {
    self == *other
  }
}

impl PartialEq<Vec<char>> for Seq {
  fn eq(&self, other: &Vec<char>) -> bool {
    self == other.as_slice()
  }
}

impl PartialEq<[u8]> for Seq {
  fn eq(&self, other: &[u8]) -> bool {
    self.data.len() == other.len() && self.data.iter().zip(other).all(|(&AsciiChar(c), &b)| c == b)
  }
}

impl PartialEq<&[u8]> for Seq {
  fn eq(&self, other: &&[u8]) -> bool {
    self == *other
  }
}

impl PartialEq<Vec<u8>> for Seq {
  fn eq(&self, other: &Vec<u8>) -> bool {
    self == other.as_slice()
  }
}

impl PartialEq<[AsciiChar]> for Seq {
  fn eq(&self, other: &[AsciiChar]) -> bool {
    self.data == other
  }
}

impl PartialEq<&[AsciiChar]> for Seq {
  fn eq(&self, other: &&[AsciiChar]) -> bool {
    self == *other
  }
}

impl PartialEq<Vec<AsciiChar>> for Seq {
  fn eq(&self, other: &Vec<AsciiChar>) -> bool {
    self == other.as_slice()
  }
}

impl core::ops::Deref for Seq {
  type Target = [AsciiChar];
  fn deref(&self) -> &Self::Target {
    &self.data
  }
}

impl core::ops::DerefMut for Seq {
  fn deref_mut(&mut self) -> &mut Self::Target {
    &mut self.data
  }
}

impl From<&Seq> for Seq {
  fn from(s: &Seq) -> Self {
    s.to_owned()
  }
}

impl From<&str> for Seq {
  fn from(s: &str) -> Self {
    Self::from_str(s)
  }
}

impl From<&[u8]> for Seq {
  fn from(slice: &[u8]) -> Self {
    Self::from_slice(slice)
  }
}

impl From<&[AsciiChar]> for Seq {
  fn from(slice: &[AsciiChar]) -> Self {
    slice.iter().copied().collect()
  }
}

impl Extend<u8> for Seq {
  fn extend<I: IntoIterator<Item = u8>>(&mut self, iter: I) {
    self.data.extend(iter.into_iter().map(AsciiChar::from));
  }
}

impl Extend<char> for Seq {
  fn extend<I: IntoIterator<Item = char>>(&mut self, iter: I) {
    self.data.extend(iter.into_iter().map(AsciiChar::from));
  }
}

impl Extend<AsciiChar> for Seq {
  fn extend<I: IntoIterator<Item = AsciiChar>>(&mut self, iter: I) {
    self.data.extend(iter);
  }
}

impl FromIterator<u8> for Seq {
  fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
    Self {
      data: iter.into_iter().map(AsciiChar::from).collect(),
    }
  }
}

impl FromIterator<AsciiChar> for Seq {
  fn from_iter<I: IntoIterator<Item = AsciiChar>>(iter: I) -> Self {
    Self {
      data: Vec::from_iter(iter),
    }
  }
}

impl FromIterator<Seq> for Seq {
  fn from_iter<I: IntoIterator<Item = Seq>>(iter: I) -> Self {
    let data = iter.into_iter().flat_map(|seq| seq.data).collect();
    Self { data }
  }
}

impl AsRef<str> for Seq {
  fn as_ref(&self) -> &str {
    self.as_str()
  }
}

impl AsRef<Seq> for Seq {
  fn as_ref(&self) -> &Seq {
    self
  }
}

impl AsRef<[AsciiChar]> for Seq {
  fn as_ref(&self) -> &[AsciiChar] {
    &self.data
  }
}

impl AsRef<[u8]> for Seq {
  fn as_ref(&self) -> &[u8] {
    self.as_str().as_bytes()
  }
}

impl core::fmt::Display for Seq {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    f.write_str(self.as_str())
  }
}

impl core::fmt::Debug for Seq {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    core::fmt::Display::fmt(self, f)
  }
}

impl core::ops::Index<usize> for Seq {
  type Output = AsciiChar;

  fn index(&self, index: usize) -> &Self::Output {
    &self.data[index]
  }
}

impl core::ops::IndexMut<usize> for Seq {
  fn index_mut(&mut self, index: usize) -> &mut Self::Output {
    &mut self.data[index]
  }
}

impl core::ops::Index<core::ops::Range<usize>> for Seq {
  type Output = [AsciiChar];

  fn index(&self, index: core::ops::Range<usize>) -> &Self::Output {
    &self.data[index]
  }
}

impl core::ops::IndexMut<core::ops::Range<usize>> for Seq {
  fn index_mut(&mut self, index: core::ops::Range<usize>) -> &mut Self::Output {
    &mut self.data[index]
  }
}

impl std::ops::Add for Seq {
  type Output = Self;

  fn add(mut self, other: Self) -> Self::Output {
    self.data.extend(other.data);
    self
  }
}

impl std::ops::Mul<usize> for Seq {
  type Output = Self;

  fn mul(mut self, rhs: usize) -> Self::Output {
    let original = self.data.clone();
    for _ in 1..rhs {
      self.data.extend(&original);
    }
    self
  }
}

impl<'a> IntoIterator for &'a Seq {
  type Item = &'a AsciiChar;
  type IntoIter = core::slice::Iter<'a, AsciiChar>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.iter()
  }
}

impl<'a> IntoIterator for &'a mut Seq {
  type Item = &'a mut AsciiChar;
  type IntoIter = core::slice::IterMut<'a, AsciiChar>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.iter_mut()
  }
}

impl IntoIterator for Seq {
  type Item = AsciiChar;
  type IntoIter = std::vec::IntoIter<AsciiChar>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.into_iter()
  }
}

#[allow(unsafe_code)]
impl std::io::Read for Seq {
  fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
    let len = std::cmp::min(buf.len(), self.len());
    // SAFETY:
    // 1. `self.data` is guaranteed to hold only ASCII characters because `Seq` enforces this invariant via `AsciiChar`.
    // 2. The length `len` is calculated as the minimum of `buf.len()` and `self.len()`, ensuring no out-of-bounds access for either slice.
    // 3. `std::ptr::copy_nonoverlapping` is safe to use here because:
    //    a. Both `self.data` and `buf` are valid, properly aligned, and non-overlapping.
    //    b. The memory regions are guaranteed to be accessible for `len` bytes.
    unsafe {
      std::ptr::copy_nonoverlapping(self.data.as_ptr().cast::<u8>(), buf.as_mut_ptr(), len);
    }
    self.data.drain(..len);
    Ok(len)
  }
}

impl std::io::Write for Seq {
  fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
    self.extend(buf.iter().copied());
    Ok(buf.len())
  }

  fn flush(&mut self) -> std::io::Result<()> {
    Ok(())
  }
}

impl serde::Serialize for Seq {
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: serde::Serializer,
  {
    serializer
      .serialize_str(self.as_str())
      .map_err(serde::ser::Error::custom)
  }
}

impl<'de> serde::Deserialize<'de> for Seq {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: serde::Deserializer<'de>,
  {
    let s = String::deserialize(deserializer)?;
    Ok(Seq::from_str(&s))
  }
}

impl JsonSchema for Seq {
  fn always_inline_schema() -> bool {
    true
  }

  fn schema_name() -> Cow<'static, str> {
    Cow::from("Seq")
  }

  fn json_schema(g: &mut SchemaGenerator) -> Schema {
    g.subschema_for::<String>()
  }
}

#[macro_export]
macro_rules! seq {
  () => (
      $crate::representation::seq::Seq::new()
  );
  ($elem:expr; $n:expr) => (
      $crate::representation::seq::Seq::from_elem($elem, $n)
  );
  ($($char:expr),* $(,)?) => {
    {
      $crate::representation::seq::Seq::from_iter([$($char),*].map(|c| c as u8).into_iter())
    }
  };
}
