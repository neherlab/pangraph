use clap::{value_parser, Args, ValueHint};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Clone, Debug, SmartDefault, Args, Serialize, Deserialize)]
pub struct AlignmentArgs {
  /// Minimum block size for alignment graph (in nucleotides)
  #[default = 100]
  #[clap(long = "len", short = 'l', default_value_t = AlignmentArgs::default().indel_len_threshold)]
  #[clap(value_hint = ValueHint::Other)]
  pub indel_len_threshold: usize,

  /// Energy cost for splitting a block during alignment merger. Controls graph fragmentation, see documentation.
  #[default = 100.0]
  #[clap(long, short = 'a', default_value_t = AlignmentArgs::default().alpha)]
  #[clap(value_hint = ValueHint::Other)]
  pub alpha: f64,

  /// Energy cost for diversity in the alignment. A high value prevents merging of distantly-related sequences in the same block, see documentation.
  #[default = 10.0]
  #[clap(long, short = 'b', default_value_t = AlignmentArgs::default().beta)]
  #[clap(value_hint = ValueHint::Other)]
  pub beta: f64,

  /// Used to set pairwise alignment sensitivity for minimap aligner. Corresponds to option -x asm5/asm10/asm20 in minimap2
  #[default = 10]
  #[clap(long, short = 's', value_parser = value_parser!(usize), default_value_t = AlignmentArgs::default().sensitivity)]
  #[clap(value_hint = ValueHint::Other)]
  pub sensitivity: usize,

  /// Sets kmer length for mmseqs2 aligner
  #[clap(long, short = 'K')]
  #[clap(value_hint = ValueHint::Other)]
  pub kmer_length: Option<usize>,

  /// Set strict maximal divergence. When toggled, proposed mergers with maximal divergence above 1/b are rejected.
  #[clap(long, short = 'S')]
  pub strict_max_divergence: bool,
}
