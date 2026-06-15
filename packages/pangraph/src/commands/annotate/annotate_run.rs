use crate::annotation::compact::{BlockCompactionStrategy, CoordinateConsensusStrategy};
use crate::annotation::lift::{LiftedAnnotation, lift_features};
use crate::annotation::matching::match_features_to_paths;
use crate::annotation::writer::{AnnotationWriter, CsvAnnotationWriter};
use crate::commands::annotate::annotate_args::{
  AnnotateCommonArgs, PangraphAnnotateArgs, PangraphAnnotateBlocksArgs, PangraphAnnotateNodesArgs,
};
use crate::io::gff::GffReader;
use crate::pangraph::pangraph::Pangraph;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;

/// Run the `annotate` command, dispatching on the chosen output granularity.
///
/// Both modes share the same front half (load graph, read GFFs, match seqids to paths, lift to
/// block-consensus coordinates); they differ only in what is written out.
pub fn annotate_run(args: PangraphAnnotateArgs) -> Result<(), Report> {
  match args {
    PangraphAnnotateArgs::Nodes(args) => annotate_run_nodes(args),
    PangraphAnnotateArgs::Blocks(args) => annotate_run_blocks(args),
  }
}

/// Load the graph, read every GFF, match seqids to paths, and lift to node-level annotations.
///
/// Returns the loaded graph alongside the lifted annotations, as block-level compaction needs the
/// graph to compute per-cluster path totals.
fn load_and_lift(common: &AnnotateCommonArgs) -> Result<(Pangraph, Vec<LiftedAnnotation>), Report> {
  let graph = Pangraph::from_path(&common.input)?;

  let mut features = Vec::new();
  for path in &common.gff {
    let read = GffReader::from_path(path)?
      .read_many()
      .wrap_err_with(|| format!("When reading GFF file: {}", path.display()))?;
    features.extend(read);
  }

  let grouped = match_features_to_paths(features, &graph, &BTreeMap::new())?;
  let lifted = lift_features(&grouped, &graph)?;

  Ok((graph, lifted))
}

/// Run `annotate nodes`: write the lossless, long-format node-level table as CSV.
fn annotate_run_nodes(args: PangraphAnnotateNodesArgs) -> Result<(), Report> {
  let PangraphAnnotateNodesArgs { common } = args;
  let (_graph, lifted) = load_and_lift(&common)?;

  let mut writer = CsvAnnotationWriter::new(&common.output, b',')?;
  writer.write_node_annotations(&lifted)?;

  Ok(())
}

/// Run `annotate blocks`: compact the node-level table into block-level consensus features (CSV).
fn annotate_run_blocks(args: PangraphAnnotateBlocksArgs) -> Result<(), Report> {
  let PangraphAnnotateBlocksArgs {
    common,
    min_frequency,
    property_threshold,
  } = args;
  let (graph, lifted) = load_and_lift(&common)?;

  let strategy = CoordinateConsensusStrategy {
    min_frequency,
    property_threshold,
  };
  let blocks = strategy.compact(&lifted, &graph)?;

  let mut writer = CsvAnnotationWriter::new(&common.output, b',')?;
  writer.write_block_annotations(&blocks)?;

  Ok(())
}
