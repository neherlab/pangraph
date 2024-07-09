use crate::align::alignment_args::AlignmentArgs;
use clap::{ArgEnum, Parser, ValueHint};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault)]
#[clap(rename = "kebab-case")]
pub enum DistanceBackend {
  #[default]
  Native,
  Mash,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault)]
#[clap(rename = "kebab-case")]
pub enum AlignmentBackend {
  #[default]
  Minimap2Lib,
  Minimap2Cli,
  Mmseqs,
}

/// Align genomes into a multiple sequence alignment graph
#[derive(Parser, Debug)]
pub struct PangraphBuildArgs {
  /// Path(s) to zero, one or multiple FASTA files with input sequences. Multiple records within one file are treated as separate genomes.
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// If no input files provided, the plain fasta input is read from standard input (stdin).
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_fastas: Vec<PathBuf>,

  #[clap(flatten, next_help_heading = "  Alignment")]
  pub aln_args: AlignmentArgs,

  /// Toggle if input genomes are circular
  #[clap(long, short = 'c')]
  pub circular: bool,

  /// Transforms all sequences to upper case
  #[clap(long, short = 'u')]
  pub upper_case: bool,

  /// Maximum number of alignment rounds to consider per pairwise graph merger
  #[clap(long, short = 'x', default_value_t = 100)]
  #[clap(value_hint = ValueHint::Other)]
  pub max_self_map: usize,

  /// Backend to use for genome similarity estimation. Similarity impacts the guide tree.
  #[clap(long, short = 'd', arg_enum, default_value_t = DistanceBackend::default())]
  #[clap(value_hint = ValueHint::Other)]
  pub distance_backend: DistanceBackend,

  /// Backend to use for pairwise genome alignment
  #[clap(long, short = 'k', arg_enum, default_value_t = AlignmentBackend::default())]
  #[clap(value_hint = ValueHint::Other)]
  pub alignment_kernel: AlignmentBackend,

  /// Random seed for block id generation
  #[clap(long)]
  pub seed: Option<u64>,
}
