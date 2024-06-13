use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Realign pancontigs of multiple sequence alignment graph
#[derive(Parser, Debug)]
pub struct PangraphPolishArgs {
  /// Path to a pangraph file (native json).
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no path provided, the uncompressed input is read from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_json: Option<PathBuf>,

  /// Maximum length: cutoff above which we won't realign
  #[clap(long, short = 'l', default_value_t = usize::MAX)]
  #[clap(value_hint = ValueHint::Other)]
  pub length: usize,

  /// Preserve case: ensure case (upper/lower) is preserved after realignment
  #[clap(long, short = 'c')]
  pub circular: bool,
}
