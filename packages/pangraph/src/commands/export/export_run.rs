use crate::commands::export::export_args::PangraphExportArgs;
use crate::io::gfa::gfa_write_file;
use crate::make_error;
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn export_run(args: PangraphExportArgs) -> Result<(), Report> {
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
    no_duplications,
    output_gfa,
    output_panx,
  } = args;

  if [&output_gfa, &output_panx].iter().all(|o| o.is_none()) {
    return make_error!("No output formats specified. Specify at least one output path.");
  }

  let graph = Pangraph::from_path(&input_json)?;

  if let Some(output_gfa) = output_gfa {
    gfa_write_file(output_gfa, &graph)?;
  }

  if let Some(output_panx) = output_panx {
    // TODO
    // panx_write(output_panx, &graph)?;
  }

  Ok(())
}
