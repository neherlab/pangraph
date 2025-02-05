use crate::commands::export::export_args::PangraphExportBlockConsensusArgs;
use crate::io::fasta::FastaWriter;
use crate::io::file::create_file_or_stdout;
use crate::pangraph::pangraph::Pangraph;
use eyre::{Report, WrapErr};

pub fn export_block_consensus(args: PangraphExportBlockConsensusArgs) -> Result<(), Report> {
  let PangraphExportBlockConsensusArgs { input_json, output } = args;
  let pangraph = Pangraph::from_path(&input_json)?;
  let mut output_fasta = FastaWriter::new(create_file_or_stdout(output)?);
  pangraph.blocks.iter().try_for_each(|(_, block)| {
    output_fasta
      .write(block.id().to_string(), &None, block.consensus())
      .wrap_err_with(|| format!("When writing consensus sequence of block {}", block.id()))
  })
}
