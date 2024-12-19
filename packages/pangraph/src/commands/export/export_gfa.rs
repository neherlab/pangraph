use crate::commands::export::export_args::PangraphExportGfaArgs;
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn export_gfa(args: PangraphExportGfaArgs) -> Result<(), Report> {
  let PangraphExportGfaArgs {
    input_json,
    output,
    minimum_length,
    maximum_length,
    minimum_depth,
    maximum_depth,
    include_sequences,
    no_duplicated,
  } = args;

  let minimum_length = minimum_length.unwrap_or(0);
  let maximum_length = maximum_length.unwrap_or(usize::MAX);
  let minimum_depth = minimum_depth.unwrap_or(0);
  let maximum_depth = maximum_depth.unwrap_or(usize::MAX);

  let pangraph_json = Pangraph::from_path(&input_json)?;

  Ok(())
}
