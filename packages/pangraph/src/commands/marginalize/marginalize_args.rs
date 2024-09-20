use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Compute all pairwise marginalizations of a multiple sequence alignment graph
#[derive(Parser, Debug)]
pub struct PangraphMarginalizeArgs {
  /// Path to multiple sequence alignment file in JSON format.
  ///
  /// Accepts plain or compressed file. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no input file provided, the uncompressed input is read from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_aln: Option<PathBuf>,

  /// Path to directory where all pairwise marginalizations will be stored
  #[clap(long, short = 'o')]
  #[clap(value_hint = ValueHint::DirPath)]
  pub output_path: PathBuf,

  /// Isolates to project onto: collapse the graph to only blocks contained by paths of the given isolates. Comma seperated list, no spaces.
  #[clap(long, short = 's', num_args=1.., use_value_delimiter = true)]
  pub strains: Vec<String>,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
