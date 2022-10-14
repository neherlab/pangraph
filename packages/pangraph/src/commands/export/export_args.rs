use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Export a pangraph to a chosen file format(s)
#[derive(Parser, Debug)]
pub struct PangraphExportArgs {
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

  /// Blocks below this length cutoff will be ignored for edges in graph
  #[clap(long, alias = "ell", default_value_t = 200)]
  #[clap(value_hint = ValueHint::Other)]
  pub edge_minimum_length: usize,

  /// Blocks below this length cutoff will be ignored for edges in graph
  #[clap(long, alias = "elu", default_value_t = usize::MAX)]
  #[clap(value_hint = ValueHint::Other)]
  pub edge_maximum_length: usize,

  /// Blocks below this depth cutoff will be ignored for edges in graph
  #[clap(long, alias = "edl", default_value_t = 0)]
  #[clap(value_hint = ValueHint::Other)]
  pub edge_minimum_depth: usize,

  /// Blocks above this depth cutoff will be ignored for edges in graph
  #[clap(long, alias = "edu", default_value_t = usize::MAX)]
  #[clap(value_hint = ValueHint::Other)]
  pub edge_maximum_depth: usize,

  /// Blocks below this length cutoff will not be exported
  #[clap(long, alias = "ll", default_value_t = 200)]
  #[clap(value_hint = ValueHint::Other)]
  pub minimum_length: usize,

  /// Blocks above this length cutoff will not be exported
  #[clap(long, alias = "lu", default_value_t = usize::MAX)]
  #[clap(value_hint = ValueHint::Other)]
  pub maximum_length: usize,

  /// Blocks below this depth cutoff will not be exported
  #[clap(long, alias = "dl", default_value_t = 0)]
  #[clap(value_hint = ValueHint::Other)]
  pub minimum_depth: usize,

  /// Blocks above this depth cutoff will not be exported
  #[clap(long, alias = "du", default_value_t = usize::MAX)]
  #[clap(value_hint = ValueHint::Other)]
  pub maximum_depth: usize,

  /// Path to directory where output will be stored
  #[clap(long, short = 'o', parse(from_os_str), default_value = "export")]
  #[clap(value_hint = ValueHint::DirPath)]
  pub output_directory: PathBuf,

  /// Basename of files
  #[clap(long, short = 'p', default_value = "pangraph")]
  #[clap(value_hint = ValueHint::DirPath)]
  pub prefix: String,

  /// Do not emit GFA file
  #[clap(long, alias = "ng")]
  pub no_export_gfa: bool,

  /// Emit vis directory to input to panX-visualization
  #[clap(long, alias = "px")]
  pub export_panx: bool,

  /// Do not export any block that contains at least one strain more than once
  #[clap(long, alias = "nd")]
  pub no_duplications: bool,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
