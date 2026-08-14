use crate::align::alignment_args::GraphMergeParams;
use clap::{Parser, ValueHint};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Merge two pangenome graphs into a single one
#[derive(Parser, Debug, SmartDefault)]
pub struct PangraphMergeArgs {
  /// Path to the first input graph, in pangraph JSON format.
  ///
  /// This graph is treated as the base: its block, node and path identifiers are preserved in the
  /// output, while those of the second graph are renumbered. When extending an existing graph with
  /// new genomes, pass the existing graph here.
  ///
  /// Accepts plain or compressed files. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`.
  /// The decompressor is chosen based on the file extension.
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub left_graph: PathBuf,

  /// Path to the second input graph, in pangraph JSON format.
  ///
  /// Its identifiers are renumbered so that they do not clash with those of the first graph. Its
  /// genomes appear after those of the first graph in the output.
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 2)]
  pub right_graph: PathBuf,

  /// Path to output JSON file with the merged pangraph.
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. If the required directory tree does not exist, it will be created.
  ///
  /// Use "-" to write the uncompressed data to standard output (stdout). This is the default, if the argument is not provided.
  #[default(PathBuf::from("-"))]
  #[clap(long, short = 'o', default_value = "-")]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub output_json: PathBuf,

  /// Sanity check: after merging verifies that every genome of the two input graphs can still be
  /// reconstructed exactly from the merged graph. Raises an error otherwise.
  #[clap(long, short = 'f')]
  pub verify: bool,

  // Declared last: it opens the "Alignment" help section, which then applies to every argument
  // registered after it, including later fields of this struct. To add one anyway, give it its own
  // `#[clap(help_heading = ...)]`, which wins. Pinned by
  // `root_args::tests::test_arguments_are_filed_under_the_expected_help_section`.
  #[clap(flatten)]
  pub merge_params: GraphMergeParams,
}
