use crate::commands::export::export_args::PangraphExportGfaArgs;
use crate::io::gfa::gfa_write_file;
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn export_gfa(args: PangraphExportGfaArgs) -> Result<(), Report> {
  let PangraphExportGfaArgs {
    input_json,
    output,
    params,
  } = args;

  let pangraph = Pangraph::from_path(&input_json)?;

  gfa_write_file(output, &pangraph, &params)?;

  Ok(())
}
