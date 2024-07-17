use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Reconstruct sequences from a pangenome graph
#[derive(Parser, Debug)]
pub struct PangraphReconstructArgs {
  /// Path to a pangenome graph file in JSON format.
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// If no input files provided, the plain fasta input is read from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_graph: Option<PathBuf>,

  /// Path to output FASTA file with reconstructed sequences.
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.
  ///
  /// Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided.
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub output_fasta: PathBuf,

  /// Path to the FASTA file with sequences to check the reconstructed sequences against. If this argument is provided, then the sequences are not being printed to standard output (stdout) as usual. Instead, if any differences are detected, a diff will be printed between the expected (original) sequence and reconstructed sequence.
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
