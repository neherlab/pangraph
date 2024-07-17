use crate::io::file::{create_file_or_stdout, open_file_or_stdin};
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use serde_json::{de::Read, Deserializer};
use std::path::Path;

/// Mitigates recursion limit error when parsing large JSONs
/// See https://github.com/serde-rs/json/issues/334
pub fn deserialize_without_recursion_limit<'de, R: Read<'de>, T: Deserialize<'de>>(
  de: &mut Deserializer<R>,
) -> Result<T, Report> {
  de.disable_recursion_limit();
  let de = serde_stacker::Deserializer::new(de);
  let obj = T::deserialize(de).wrap_err("When parsing JSON")?;
  Ok(obj)
}

pub fn json_read_str<T: for<'de> Deserialize<'de>>(s: &str) -> Result<T, Report> {
  let mut de = Deserializer::from_str(s);
  deserialize_without_recursion_limit(&mut de)
}

pub fn json_read_bytes<T: for<'de> Deserialize<'de>>(bytes: &[u8]) -> Result<T, Report> {
  let mut de = Deserializer::from_slice(bytes);
  deserialize_without_recursion_limit(&mut de)
}

pub fn json_read_file<T: for<'de> Deserialize<'de>, P: AsRef<Path>>(filepath: &Option<P>) -> Result<T, Report> {
  let r = open_file_or_stdin(filepath)?;
  json_read_reader(r)
}

pub fn json_read_reader<T: for<'de> Deserialize<'de>>(r: impl std::io::Read) -> Result<T, Report> {
  let mut de = Deserializer::from_reader(r);
  deserialize_without_recursion_limit(&mut de)
}

pub fn json_stringify<T: Serialize>(obj: &T) -> Result<String, Report> {
  serde_json::to_string_pretty(obj).wrap_err("When converting an entry to JSON string")
}

pub fn json_write<T: Serialize>(filepath: impl AsRef<Path>, obj: &T) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let file = create_file_or_stdout(filepath)?;
  serde_json::to_writer_pretty(file, obj).wrap_err("When writing JSON to file: {filepath:#?}")
}
