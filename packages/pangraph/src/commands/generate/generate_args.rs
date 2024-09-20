use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Output a simulated sequence alignment graph
#[derive(Parser, Debug)]
pub struct PangraphGenerateArgs {
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

  /// Mutation rate: rate of mutations per site per genome per generation
  #[clap(long, short = 'm', default_value_t = 1.0e-5)]
  #[clap(value_hint = ValueHint::Other)]
  pub snp_rate: f64,

  /// HGT rate: rate of horizontal gene transfer events per genome per generation
  #[clap(long, short = 'r', default_value_t = 0.0)]
  #[clap(value_hint = ValueHint::Other)]
  pub hgt_rate: f64,

  /// Length variance prefactor (σ=L/s): divisor of mean length that sets the variance
  #[clap(long, short = 's', default_value_t = 10.0)]
  #[clap(value_hint = ValueHint::Other)]
  pub sigma_pre: f64,

  /// Deletion rate: Rate of deletion events per genome per generation
  #[clap(long, short = 'd', default_value_t = 0.0)]
  #[clap(value_hint = ValueHint::Other)]
  pub delete_rate: f64,

  /// Inversion rate: rate of inversion events per genome per generation
  #[clap(long, short = 'i', default_value_t = 0.0)]
  #[clap(value_hint = ValueHint::Other)]
  pub invert_rate: f64,

  /// Path where to output graph
  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::FilePath)]
  pub graph_output: PathBuf,

  /// Generations to simulate: number of generations simulated under WF model
  #[clap(long, short = 't', default_value_t = 35)]
  #[clap(value_hint = ValueHint::Other)]
  pub time: usize,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
