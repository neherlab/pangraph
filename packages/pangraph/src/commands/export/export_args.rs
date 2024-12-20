use clap::{Parser, Subcommand, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Export a pangraph to a chosen file format
#[derive(Subcommand, Debug, Clone)]
#[clap(verbatim_doc_comment)]
pub enum PangraphExportArgs {
  /// Export to GFA v1 format
  ///
  /// See GFA v1 format specifications: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
  Gfa(PangraphExportGfaArgs),

  /// Export block consensus sequences to a fasta file
  BlockConsensus(PangraphExportBlockConsensusArgs),

  /// Export aligned or unaligned sequences for each block. Note that alignments exclude insertions
  BlockSequences(PangraphExportBlockSequencesArgs),

  /// Export the core-genome alignment
  CoreGenome(PangraphExportCoreAlignmentArgs),
}

#[derive(Parser, Debug, Clone)]
pub struct PangraphExportGfaArgs {
  /// Path to a pangraph file (native json).
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no path provided, the uncompressed input is read from standard input (stdin).
  #[clap(display_order = 1, value_hint = ValueHint::FilePath)]
  pub input_json: Option<PathBuf>,

  /// Path to output GFA file.
  ///
  /// See GFA v1 format specifications: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.
  ///
  /// Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided.
  #[clap(long, short = 'o', default_value = "-", value_hint = ValueHint::AnyPath)]
  pub output: PathBuf,

  /// Blocks below this length cutoff will not be exported
  #[clap(long, value_hint = ValueHint::Other)]
  pub minimum_length: Option<usize>,

  /// Blocks above this length cutoff will not be exported
  #[clap(long, value_hint = ValueHint::Other)]
  pub maximum_length: Option<usize>,

  /// Blocks below this depth cutoff will not be exported
  #[clap(long, value_hint = ValueHint::Other)]
  pub minimum_depth: Option<usize>,

  /// Blocks above this depth cutoff will not be exported
  #[clap(long, value_hint = ValueHint::Other)]
  pub maximum_depth: Option<usize>,

  /// Include block sequences in the GFA file
  #[clap(long)]
  pub include_sequences: bool,

  /// Exclude blocks that are duplicated in any path
  #[clap(long)]
  pub no_duplicated: bool,
}

#[derive(Parser, Debug, Clone)]
pub struct PangraphExportBlockConsensusArgs {
  /// Path to a pangraph file (native json).
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no path provided, the uncompressed input is read from standard input (stdin).
  #[clap(display_order = 1, value_hint = ValueHint::FilePath)]
  pub input_json: Option<PathBuf>,

  /// Path to output FASTA file.
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.
  ///
  /// Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided.
  #[clap(long, short = 'o', default_value = "-", value_hint = ValueHint::AnyPath)]
  pub output: PathBuf,
}

#[derive(Parser, Debug, Clone)]
pub struct PangraphExportBlockSequencesArgs {
  /// Path to a pangraph file (native json).
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no path provided, the uncompressed input is read from standard input (stdin).
  #[clap(display_order = 1, value_hint = ValueHint::FilePath)]
  pub input_json: Option<PathBuf>,

  /// Path to directory to write output FASTA files to
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[clap(long, short = 'o', value_hint = ValueHint::AnyPath)]
  pub output: PathBuf,

  /// If set, then the full block sequences are exported but not aligned.
  #[clap(long)]
  pub unaligned: bool,
}

#[derive(Parser, Debug, Clone)]
pub struct PangraphExportCoreAlignmentArgs {
  /// Path to a pangraph file (native json).
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no path provided, the uncompressed input is read from standard input (stdin).
  #[clap(display_order = 1, value_hint = ValueHint::FilePath)]
  pub input_json: Option<PathBuf>,

  /// Path to output FASTA file.
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.
  ///
  /// Use "-" to write the uncompressed to standard output (stdout). This is the default, if the argument is not provided.
  #[clap(long, short = 'o', default_value = "-", value_hint = ValueHint::AnyPath)]
  pub output: PathBuf,

  /// Specify the strain to use as a reference for the alignment. Core blocks are ordered and oriented (forward or reverse) according to the reference strain.
  #[clap(long, value_hint = ValueHint::Other)]
  pub guide_strain: Option<String>,

  /// If set, then the full core sequences are exported but not aligned.
  ///
  /// They should be linearly alignable and can be fed to an external aligner.
  #[clap(long)]
  pub unaligned: bool,
}
