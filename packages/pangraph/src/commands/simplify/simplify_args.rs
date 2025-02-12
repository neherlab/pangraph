use clap::{ArgAction, Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Generates a simplified graph that only contains a subset of the input genomes.
#[derive(Parser, Debug)]
pub struct PangraphSimplifyArgs {
  /// Path to Pangraph JSON.
  ///
  /// Accepts plain or compressed file. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no input file provided, the uncompressed input is read from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input: Option<PathBuf>,

  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub output: PathBuf,

  /// Isolates to project onto: collapse the graph to only blocks contained by paths of the given isolates. List of strain names, comma-delimited without spaces.
  #[clap(long, short = 's', required = true, num_args = 1, action = ArgAction::Set, value_delimiter=',')]
  // See: https://github.com/clap-rs/clap/issues/4942#issuecomment-1565139247
  pub strains: Vec<String>,
}
