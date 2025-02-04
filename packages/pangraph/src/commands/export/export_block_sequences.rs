use crate::commands::export::export_args::PangraphExportBlockSequencesArgs;
use crate::io::fasta::FastaWriter;
use crate::io::file::create_file_or_stdout;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::RecordNaming;
use clap::Parser;
use eyre::{Context, Report};

#[derive(Parser, Debug, Default, Clone)]
pub struct ExportBlockSequencesParams {
  /// If set, then the full non-aligned block sequences are exported.
  #[clap(long)]
  pub unaligned: bool,
}

pub fn export_block_sequences(args: PangraphExportBlockSequencesArgs) -> Result<(), Report> {
  let PangraphExportBlockSequencesArgs {
    input_json,
    output,
    params,
  } = args;

  let pangraph = Pangraph::from_path(&input_json)?;

  pangraph.blocks.iter().try_for_each(|(_, block)| {
    let fasta_filepath = output.join(format!("block_{}.fa", block.id()));

    let mut output_fasta = FastaWriter::new(create_file_or_stdout(&fasta_filepath)?);

    block
      .sequences(&pangraph, !params.unaligned, RecordNaming::Node)
      .enumerate()
      .try_for_each(|(index, (id, seq))| {
        {
          let seq = seq?;
          output_fasta.write(&id, &None, &seq)
        }
        .wrap_err_with(|| format!("When writing sequence #{index} '{id}'"))
        .wrap_err_with(|| format!("When writing sequences of block {}", block.id()))
      })
  })
}
