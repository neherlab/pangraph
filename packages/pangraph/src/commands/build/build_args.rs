use crate::pangraph::split_matches::SplitMatchesArgs;
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
  Minimap2,
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

  #[clap(flatten)]
  pub split_matches_args: SplitMatchesArgs,

  /// Energy cost for introducing junction due to alignment merger
  #[clap(long, short = 'a', default_value_t = 100.0)]
  #[clap(value_hint = ValueHint::Other)]
  pub alpha: f64,

  /// Energy cost for interblock diversity due to alignment merger
  #[clap(long, short = 'b', default_value_t = 10.0)]
  #[clap(value_hint = ValueHint::Other)]
  pub beta: f64,

  /// Toggle if input genomes are circular
  #[clap(long, short = 'c')]
  pub circular: bool,

  /// Transforms all sequences to upper case
  #[clap(long, short = 'u')]
  pub upper_case: bool,

  /// Used to set pairwise alignment sensitivity
  #[clap(long, short = 's', possible_values(&["5", "10", "20"]), default_value_t = 10)]
  #[clap(value_hint = ValueHint::Other)]
  pub sensitivity: usize,

  /// Maximum number of self mappings to consider per pairwise graph merger
  #[clap(long, short = 'x', default_value_t = 100)]
  #[clap(value_hint = ValueHint::Other)]
  pub max_self_map: usize,

  /// Backend to use for pairwise genome alignment
  #[clap(long, short = 'd', arg_enum, default_value_t = DistanceBackend::default())]
  #[clap(value_hint = ValueHint::Other)]
  pub distance_backend: DistanceBackend,

  /// Backend to use for pairwise genome alignment
  #[clap(long, short = 'k', arg_enum, default_value_t = AlignmentBackend::default())]
  #[clap(value_hint = ValueHint::Other)]
  pub alignment_kernel: AlignmentBackend,

  #[clap(long, short = 'K')]
  #[clap(value_hint = ValueHint::Other)]
  pub kmer_length: Option<usize>,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
