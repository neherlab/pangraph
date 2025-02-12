use crate::io::json::json_or_yaml_write_file;
use eyre::Report;
use schemars::{schema_for, JsonSchema};
use std::path::Path;

/// Write JSON Schema derived from type T into a file
pub fn jsonschema_write_file<T: JsonSchema>(output_file: &Option<impl AsRef<Path>>) -> Result<(), Report> {
  let schema = schema_for!(T);
  let output_file = output_file.as_ref().map_or_else(|| Path::new("-"), AsRef::as_ref);
  json_or_yaml_write_file(output_file, &schema)
}
