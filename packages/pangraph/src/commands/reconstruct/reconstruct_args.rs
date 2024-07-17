use clap::{ArgEnum, Parser, ValueHint};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// reconstruct sequences from a pangenome graph
#[derive(Parser, Debug)]
pub struct PangraphReconstructArgs {
  /// Path to a pangenome graph file in JSON format.
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_graph: PathBuf,
}
