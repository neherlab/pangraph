use crate::commands::polish::polish_args::PangraphPolishArgs;
use crate::io::json::json_write;
use crate::io::pangraph_json::PangraphJson;
use crate::utils::random::get_random_number_generator;
use eyre::Report;

pub fn polish_run(args: &PangraphPolishArgs) -> Result<(), Report> {
  let PangraphPolishArgs {
    input_json,
    length,
    circular,
    seed,
  } = &args;

  let rng = get_random_number_generator(seed);

  let input_pangraph_json = PangraphJson::from_path(input_json)?;
  let output_pangraph_json = polish(&input_pangraph_json, args)?;
  json_write("-", &output_pangraph_json)?;

  Ok(())
}

pub fn polish(input_pangraph_json: &PangraphJson, args: &PangraphPolishArgs) -> Result<PangraphJson, Report> {
  // TODO: create proper output pangraph JSON ;)
  let output_pangraph_json = input_pangraph_json.clone();
  Ok(output_pangraph_json)
}
