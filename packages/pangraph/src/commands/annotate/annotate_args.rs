use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Lift genome annotations onto the pangenome graph.
///
/// Reads one or more GFF3 annotation files and places each feature on the graph node(s) it overlaps,
/// translating its coordinates into block-consensus coordinates. The result is a long-format,
/// node-level table (the lossless source of truth), written as CSV.
///
/// Annotation `seqid`s are matched to graph path names by exact string equality; any `seqid` that
/// does not correspond to a path is a hard error (annotation seqids must match the FASTA record
/// names used to build the graph).
#[derive(Parser, Debug)]
pub struct PangraphAnnotateArgs {
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

  /// Path to a GFF3 annotation file. Repeat the flag to provide multiple files.
  ///
  /// Accepts plain or compressed files (`gz`, `bz2`, `xz`, `zstd`), chosen by file extension. At
  /// least one file is required. Annotation `seqid`s must match the graph path names exactly.
  #[clap(long = "gff", required = true, value_hint = ValueHint::FilePath)]
  #[clap(display_order = 2)]
  pub gff: Vec<PathBuf>,

  /// Path to the output node-level annotation table (CSV).
  ///
  /// Will be created if it does not exist. The output is compressed if the path ends in a known
  /// compression extension (`gz`, `bz2`, `xz`, `zstd`). Use `-` to write uncompressed CSV to
  /// standard output (stdout).
  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub output: PathBuf,
}
