use eyre::Report;
use pangraph::commands::schema::generate_schema::{generate_schema, PangraphGenerateSchemaArgs};
use std::path::PathBuf;

fn main() -> Result<(), Report> {
  generate_schema(&PangraphGenerateSchemaArgs {
    output_pangraph_schema: PathBuf::from("../pangraph-schemas/Pangraph.schema.json"),
  })?;
  Ok(())
}