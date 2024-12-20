use crate::commands::export::export_args::PangraphExportBlockSequencesArgs;
use crate::io::fasta::FastaWriter;
use crate::io::file::create_file_or_stdout;
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn export_block_sequences(args: PangraphExportBlockSequencesArgs) -> Result<(), Report> {
  let PangraphExportBlockSequencesArgs {
    input_json,
    output,
    unaligned,
  } = args;

  let pangraph_json = Pangraph::from_path(&input_json)?;
  let output_fasta = FastaWriter::new(create_file_or_stdout(output)?);

  Ok(())
}
