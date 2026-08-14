use crate::align::alignment_args::GraphMergeParams;
use clap::{Parser, ValueHint};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Align genomes into a pangenome graph
#[derive(Parser, Debug, SmartDefault)]
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

  /// Path to output JSON file with resulting pangraph.
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.
  ///
  /// Use "-" to write the uncompressed data to standard output (stdout). This is the default, if the argument is not provided.
  #[default(PathBuf::from("-"))]
  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub output_json: PathBuf,

  /// Toggle if input genomes are circular
  #[clap(long, short = 'c')]
  pub circular: bool,

  /// Sanity check: after construction verifies that the original sequences can be reconstructed exactly from the resulting pangraph. Raises an error otherwise.
  #[clap(long, short = 'f')]
  pub verify: bool,

  /// Toggle to disable progress bar. Notice that the progress bar is only displayed if the output is specified via the `-o` argument.
  #[clap(long)]
  pub no_progress_bar: bool,

  /// Path to a Newick-format guide tree to use instead of the default neighbor-joining tree.
  ///
  /// When provided, the tree's topology drives the bottom-up graph-merging order. Each input
  /// FASTA sequence must appear exactly once as a leaf (matched by sequence name), and every
  /// internal node must be strictly bifurcating. Branch lengths and internal labels, if present,
  /// are ignored. Accepts plain or compressed files (gz, bz2, xz, zst).
  #[clap(long, value_hint = ValueHint::FilePath)]
  pub guide_tree: Option<PathBuf>,

  // Declared last: it opens the "Alignment" help section, which then applies to every argument
  // registered after it, including later fields of this struct. To add one anyway, give it its own
  // `#[clap(help_heading = ...)]`, which wins. Pinned by
  // `root_args::tests::test_arguments_are_filed_under_the_expected_help_section`.
  #[clap(flatten)]
  pub merge_params: GraphMergeParams,
}
