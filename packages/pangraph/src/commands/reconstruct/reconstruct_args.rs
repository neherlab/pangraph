use clap::{Parser, ValueHint};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Reconstruct sequences from a pangenome graph
#[derive(Parser, Debug, SmartDefault)]
pub struct PangraphReconstructArgs {
  /// Path to a pangenome graph file in JSON format.
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no input file is provided, the plain JSON input is read from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_graph: Option<PathBuf>,

  /// Path to output FASTA file with reconstructed sequences.
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.
  ///
  /// Use "-" to write the uncompressed data to standard output (stdout). This is the default, if the argument is not provided.
  ///
  /// Records are written in order of path id, which reproduces the order of the original input FASTA
  /// only for graphs produced directly by `pangraph build`. A graph produced by `pangraph merge`
  /// renumbers its path ids, so consumers should match records by genome name rather than by
  /// position.
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub output_fasta: PathBuf,

  /// Path to the FASTA file with sequences to check the reconstructed sequences against. If this argument is provided, then the sequences are not written out as usual: nothing is produced on success, and the first difference found is reported as an error.
  ///
  /// Genomes are matched by name, so the order of the records is irrelevant. The file must contain
  /// exactly the genomes of the graph: a record the graph does not contain, a genome missing from
  /// the file, or a repeated name are all errors. Every path of the graph must be named.
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// Use "-" to read uncompressed FASTA from standard input (stdin).
  #[clap(long, short = 'f')]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub verify: Option<PathBuf>,
}
