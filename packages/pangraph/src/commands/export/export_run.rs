use crate::commands::export::export_args::PangraphExportArgs;
use crate::commands::export::export_block_consensus::export_block_consensus;
use crate::commands::export::export_block_sequences::export_block_sequences;
use crate::commands::export::export_core_genome::export_core_genome;
use crate::commands::export::export_gfa::export_gfa;
use eyre::Report;

pub fn export_run(args: PangraphExportArgs) -> Result<(), Report> {
  match args {
    PangraphExportArgs::Gfa(args) => {
      export_gfa(args)?;
    },
    PangraphExportArgs::BlockConsensus(args) => {
      export_block_consensus(args)?;
    },
    PangraphExportArgs::BlockSequences(args) => {
      export_block_sequences(args)?;
    },
    PangraphExportArgs::CoreGenome(args) => {
      export_core_genome(args)?;
    },
  }
  Ok(())
}
