use crate::commands::polish::polish_args::PangraphPolishArgs;
use crate::io::json::{json_write_file, JsonPretty};
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn polish_run(args: &PangraphPolishArgs) -> Result<(), Report> {
  let PangraphPolishArgs {
    input_json,
    length,
    circular,
  } = &args;

  let input_pangraph_json = Pangraph::from_path(input_json)?;
  let output_pangraph_json = polish(&input_pangraph_json, args)?;
  json_write_file("-", &output_pangraph_json, JsonPretty(true))?;

  Ok(())
}

pub fn polish(input_pangraph_json: &Pangraph, args: &PangraphPolishArgs) -> Result<Pangraph, Report> {
  // TODO: create proper output pangraph JSON ;)
  let output_pangraph_json = input_pangraph_json.clone();
  Ok(output_pangraph_json)
}
