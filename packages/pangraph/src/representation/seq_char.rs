use schemars::{JsonSchema, Schema, SchemaGenerator};
use std::borrow::Cow;
use std::fmt::Write as StdFmtWrite;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AsciiChar(pub u8);

impl AsciiChar {
  pub const fn new(value: u8) -> Self {
    Self(value)
  }

  pub const fn inner(&self) -> u8 {
    self.0
  }

  pub fn from_str(s: &str) -> Self {
    debug_assert!(s.is_ascii());
    debug_assert!(s.len() == 1);
    Self(s.as_bytes()[0])
  }
}

impl core::fmt::Display for AsciiChar {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    f.write_char(self.0 as char)
  }
}

impl core::fmt::Debug for AsciiChar {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    core::fmt::Display::fmt(self, f)
  }
}

impl From<u8> for AsciiChar {
  fn from(item: u8) -> Self {
    AsciiChar(item)
  }
}

impl From<u16> for AsciiChar {
  fn from(item: u16) -> Self {
    AsciiChar(item as u8)
  }
}

impl From<u32> for AsciiChar {
  fn from(item: u32) -> Self {
    AsciiChar(item as u8)
  }
}

impl From<u64> for AsciiChar {
  fn from(item: u64) -> Self {
    AsciiChar(item as u8)
  }
}

impl From<usize> for AsciiChar {
  fn from(item: usize) -> Self {
    AsciiChar(item as u8)
  }
}

impl From<char> for AsciiChar {
  fn from(item: char) -> Self {
    AsciiChar(item as u8)
  }
}

impl From<AsciiChar> for u8 {
  fn from(item: AsciiChar) -> Self {
    item.0
  }
}

impl From<AsciiChar> for u16 {
  fn from(item: AsciiChar) -> Self {
    item.0 as u16
  }
}

impl From<AsciiChar> for u32 {
  fn from(item: AsciiChar) -> Self {
    item.0 as u32
  }
}

impl From<AsciiChar> for u64 {
  fn from(item: AsciiChar) -> Self {
    item.0 as u64
  }
}

impl From<AsciiChar> for usize {
  fn from(item: AsciiChar) -> Self {
    item.0 as usize
  }
}

impl From<AsciiChar> for char {
  fn from(item: AsciiChar) -> Self {
    item.0 as char
  }
}

impl serde::Serialize for AsciiChar {
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: serde::Serializer,
  {
    serializer
      .serialize_str(&self.to_string())
      .map_err(serde::ser::Error::custom)
  }
}

impl<'de> serde::Deserialize<'de> for AsciiChar {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: serde::Deserializer<'de>,
  {
    let s = String::deserialize(deserializer)?;
    Ok(AsciiChar::from_str(&s))
  }
}

impl JsonSchema for AsciiChar {
  fn always_inline_schema() -> bool {
    true
  }

  fn schema_name() -> Cow<'static, str> {
    Cow::from("AsciiChar")
  }

  fn json_schema(g: &mut SchemaGenerator) -> Schema {
    g.subschema_for::<char>()
  }
}
