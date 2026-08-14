use clap::{Args, ValueEnum, ValueHint, value_parser};
use color_eyre::owo_colors::{AnsiColors, OwoColorize};
use color_eyre::{Help, SectionExt};
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use strum_macros::Display;

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
}

#[derive(
  Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, SmartDefault, Display, Serialize, Deserialize,
)]
#[clap(rename_all = "kebab-case")]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum AlignmentBackend {
  #[default]
  Minimap2,
  Mmseqs,
}

/// Parameters that govern how two graphs are merged into one.
///
/// Shared by `pangraph build` (which merges at every node of the guide tree) and
/// `pangraph merge` (which performs a single such merger), so that both commands expose an
/// identical set of alignment options.
#[derive(Clone, Debug, SmartDefault, Args, Serialize, Deserialize)]
pub struct GraphMergeParams {
  #[clap(flatten, next_help_heading = "Alignment")]
  pub aln_args: AlignmentArgs,

  /// Maximum number of alignment rounds to consider per pairwise graph merger
  #[default = 100]
  #[clap(long, short = 'x', default_value_t = GraphMergeParams::default().max_self_map)]
  #[clap(value_hint = ValueHint::Other)]
  pub max_self_map: usize,

  /// Backend to use for pairwise genome alignment
  ///
  /// Nb: `mmseqs` is more sensitive to highly-diverged sequences, but slower and requires more memory.
  /// It is not provided with Pangraph, so you need to install it separately (see: https://github.com/soedinglab/MMseqs2)
  #[clap(long, short = 'k',  default_value_t = GraphMergeParams::default().alignment_kernel)]
  #[clap(value_hint = ValueHint::Other)]
  pub alignment_kernel: AlignmentBackend,

  /// For within-block alignment: excess bandwidth for internal stripes.
  /// Can be increased to improve block alignment quality, at the cost of computation time and memory usage.
  #[default = 5]
  #[clap(long, default_value_t = GraphMergeParams::default().extra_band_width)]
  #[clap(value_hint = ValueHint::Other)]
  pub extra_band_width: usize,

  /// For within-block alignment: number of times Nextclade will retry alignment with more relaxed results if alignment band boundaries are hit.
  #[default = 4]
  #[clap(long, default_value_t = GraphMergeParams::default().max_alignment_attempts)]
  #[clap(value_hint = ValueHint::Other)]
  pub max_alignment_attempts: usize,
}

/// Checks that the selected alignment backend is usable before any real work is started.
/// Only the external backends need checking: the minimap2 one is built into pangraph.
pub fn check_alignment_backend_available(params: &GraphMergeParams) -> Result<(), Report> {
  if params.alignment_kernel == AlignmentBackend::Mmseqs {
    // check that mmseqs is available in PATH
    std::process::Command::new("mmseqs")
      .arg("--help")
      .output()
      .wrap_err("When executing `mmseqs --help`")
      .section(
        "Please make sure that `mmseqs` is installed, available in PATH and is functional outside of pangraph. For more details, refer to mmseqs documentation at https://github.com/soedinglab/MMseqs2"
          .color(AnsiColors::Cyan)
          .header("Suggestion:"),
      )?;
  }

  Ok(())
}
