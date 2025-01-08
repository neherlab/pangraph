use crate::commands::export::export_block_sequences::ExportBlockSequencesParams;
use crate::commands::export::export_core_genome::ExportCoreAlignmentParams;
use crate::io::gfa::GfaWriteParams;
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

  /// Export aligned or unaligned sequences for each block in separate fasta files. Note that alignments exclude insertions.
  BlockSequences(PangraphExportBlockSequencesArgs),

  /// Export the core-genome alignment. Note that alignment excludes insertions.
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

  #[clap(flatten)]
  pub params: GfaWriteParams,
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

  #[clap(flatten)]
  pub params: ExportBlockSequencesParams,
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

  #[clap(flatten)]
  pub params: ExportCoreAlignmentParams,
}
