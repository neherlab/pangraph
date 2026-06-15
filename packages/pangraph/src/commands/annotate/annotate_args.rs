use clap::{Parser, Subcommand, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

/// Lift genome annotations onto the pangenome graph.
///
/// Reads one or more GFF3 annotation files and places each feature on the graph node(s) it overlaps,
/// translating its coordinates into block-consensus coordinates. Choose the granularity of the
/// output with a subcommand:
///
/// - `nodes`  — the lossless, long-format node-level table (one row per feature/overlapped node).
/// - `blocks` — block-level consensus features, collapsing the redundant per-node placements that
///   most genomes share into one row per cluster.
///
/// Annotation `seqid`s are matched to graph path names by exact string equality; any `seqid` that
/// does not correspond to a path is a hard error (annotation seqids must match the FASTA record
/// names used to build the graph).
#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum PangraphAnnotateArgs {
  /// Lift annotations to per-node block-consensus coordinates (lossless, long-format CSV).
  Nodes(PangraphAnnotateNodesArgs),

  /// Compact node-level annotations into block-level consensus features (CSV).
  Blocks(PangraphAnnotateBlocksArgs),
}

/// Options shared by every `annotate` subcommand: the graph, the GFF inputs, and the output path.
#[derive(Parser, Debug)]
pub struct AnnotateCommonArgs {
  /// Path to Pangraph JSON.
  ///
  /// Accepts plain or compressed file. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension.
  ///
  /// If no input file provided, the uncompressed input is read from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input: Option<PathBuf>,

  /// Path to a GFF3 annotation file. Repeat the flag to provide multiple files.
  ///
  /// Accepts plain or compressed files (`gz`, `bz2`, `xz`, `zstd`), chosen by file extension. At
  /// least one file is required. Annotation `seqid`s must match the graph path names exactly.
  #[clap(long = "gff", required = true, value_hint = ValueHint::FilePath)]
  #[clap(display_order = 2)]
  pub gff: Vec<PathBuf>,

  /// Path to the output annotation table (CSV).
  ///
  /// Will be created if it does not exist. The output is compressed if the path ends in a known
  /// compression extension (`gz`, `bz2`, `xz`, `zstd`). Use `-` to write uncompressed CSV to
  /// standard output (stdout).
  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub output: PathBuf,
}

/// Arguments for `annotate nodes`: produce the lossless node-level table.
#[derive(Parser, Debug)]
pub struct PangraphAnnotateNodesArgs {
  #[clap(flatten)]
  pub common: AnnotateCommonArgs,
}

/// Arguments for `annotate blocks`: compact the node-level table into block-level consensus features.
///
/// The two thresholds mirror the fields of `CoordinateConsensusStrategy`; their defaults are kept in
/// sync with that type's `Default` (0.9 / 0.5).
#[derive(Parser, Debug)]
pub struct PangraphAnnotateBlocksArgs {
  #[clap(flatten)]
  pub common: AnnotateCommonArgs,

  /// Minimum frequency required to emit a block-level cluster.
  ///
  /// A cluster is kept when the number of supporting genomes `M >= ceil(min_frequency * N)`, where
  /// `N` is the number of paths traversing the cluster's block(s) (so a gene is not penalised for
  /// being absent in genomes that lack the block entirely).
  #[clap(long, default_value_t = 0.9, value_parser = parse_fraction)]
  #[clap(value_hint = ValueHint::Other)]
  pub min_frequency: f64,

  /// Minimum supporter agreement required to promote a consensus name or attribute value.
  ///
  /// For each cluster, a `name`/attribute value is written only if at least this fraction of the
  /// supporting genomes agree on it; otherwise the field is left empty.
  #[clap(long, default_value_t = 0.5, value_parser = parse_fraction)]
  #[clap(value_hint = ValueHint::Other)]
  pub property_threshold: f64,
}

/// Parse a threshold given as a fraction in the closed unit interval `[0, 1]`.
///
/// Both `annotate blocks` thresholds are fractions; a value outside `[0, 1]` is always a mistake (it
/// would silently emit nothing, or promote every value), so it is rejected at parse time rather than
/// failing quietly downstream. `NaN` and infinities fall outside the range and are rejected too.
fn parse_fraction(s: &str) -> Result<f64, String> {
  let value: f64 = s.parse().map_err(|err| format!("`{s}` is not a valid number: {err}"))?;
  if (0.0..=1.0).contains(&value) {
    Ok(value)
  } else {
    Err(format!(
      "must be a fraction between 0 and 1 (inclusive), but got `{value}`"
    ))
  }
}

#[cfg(test)]
mod tests {
  use super::parse_fraction;

  #[test]
  fn parse_fraction_accepts_closed_unit_interval() {
    parse_fraction("0").unwrap();
    parse_fraction("0.5").unwrap();
    parse_fraction("1").unwrap();
  }

  #[test]
  fn parse_fraction_rejects_out_of_range_and_non_numeric() {
    parse_fraction("-0.01").unwrap_err();
    parse_fraction("1.01").unwrap_err();
    parse_fraction("NaN").unwrap_err();
    parse_fraction("inf").unwrap_err();
    parse_fraction("abc").unwrap_err();
  }
}
