use clap::{Args, ValueHint};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Clone, Debug, SmartDefault, Args, Serialize, Deserialize)]
pub struct AlignmentArgs {
  /// Minimum block size for alignment graph (in nucleotides)
  #[default = 100]
  #[clap(long = "len", short = 'l', default_value_t = AlignmentArgs::default().indel_len_threshold)]
  #[clap(value_hint = ValueHint::Other)]
  pub indel_len_threshold: usize,

  /// Energy cost for introducing junction due to alignment merger
  #[default = 100.0]
  #[clap(long, short = 'a', default_value_t = AlignmentArgs::default().alpha)]
  #[clap(value_hint = ValueHint::Other)]
  pub alpha: f64,

  /// Energy cost for interblock diversity due to alignment merger
  #[default = 10.0]
  #[clap(long, short = 'b', default_value_t = AlignmentArgs::default().beta)]
  #[clap(value_hint = ValueHint::Other)]
  pub beta: f64,

  /// Used to set pairwise alignment sensitivity for minimap aligner. Corresponds to option -x asm5/asm10/asm20 in minimap2
  #[default = 10]
  #[clap(long, short = 's', possible_values(&["5", "10", "20"]), default_value_t = AlignmentArgs::default().sensitivity)]
  #[clap(value_hint = ValueHint::Other)]
  pub sensitivity: usize,

  /// Sets kmer length for mmseqs2 aligner
  #[clap(long, short = 'K')]
  #[clap(value_hint = ValueHint::Other)]
  pub kmer_length: Option<usize>,
}
