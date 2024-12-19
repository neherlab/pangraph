use crate::commands::export::export_args::PangraphExportCoreAlignmentArgs;
use crate::io::fasta::FastaWriter;
use crate::io::file::create_file_or_stdout;
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn export_core_genome(args: PangraphExportCoreAlignmentArgs) -> Result<(), Report> {
  let PangraphExportCoreAlignmentArgs {
    input_json,
    output,
    guide_strain,
    unaligned,
  } = args;

  let pangraph_json = Pangraph::from_path(&input_json)?;
  let output_fasta = FastaWriter::new(create_file_or_stdout(output)?);

  Ok(())
}
