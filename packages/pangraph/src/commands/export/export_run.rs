use crate::commands::export::export_args::PangraphExportArgs;
use crate::pangraph::pangraph::Pangraph;
use crate::utils::random::get_random_number_generator;
use eyre::Report;

pub fn export_run(args: &PangraphExportArgs) -> Result<(), Report> {
  let PangraphExportArgs {
    input_json,
    edge_minimum_length,
    edge_maximum_length,
    edge_minimum_depth,
    edge_maximum_depth,
    minimum_length,
    maximum_length,
    minimum_depth,
    maximum_depth,
    output_directory,
    prefix,
    no_export_gfa,
    export_panx,
    no_duplications,
    seed,
  } = &args;

  let rng = get_random_number_generator(seed);

  let pangraph_json = Pangraph::from_path(input_json)?;

  Ok(())
}
