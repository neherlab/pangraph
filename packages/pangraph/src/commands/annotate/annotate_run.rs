use crate::annotation::lift::lift_features;
use crate::annotation::matching::match_features_to_paths;
use crate::annotation::writer::{AnnotationWriter, CsvAnnotationWriter};
use crate::commands::annotate::annotate_args::PangraphAnnotateArgs;
use crate::io::gff::GffReader;
use crate::pangraph::pangraph::Pangraph;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;

/// Run the `annotate` command: lift GFF features onto graph nodes and write the node-level table.
///
/// Loads the graph, reads every GFF file into the internal `Feature` model, matches each feature's
/// `seqid` to a graph path (exact match; unmatched seqids are a hard error), lifts the features to
/// block-consensus coordinates, and writes the resulting node-level annotations as CSV.
pub fn annotate_run(args: PangraphAnnotateArgs) -> Result<(), Report> {
  let PangraphAnnotateArgs { input, gff, output } = args;

  let graph = Pangraph::from_path(&input)?;

  let mut features = Vec::new();
  for path in &gff {
    let read = GffReader::from_path(path)?
      .read_many()
      .wrap_err_with(|| format!("When reading GFF file: {}", path.display()))?;
    features.extend(read);
  }

  let grouped = match_features_to_paths(features, &graph, &BTreeMap::new())?;
  let lifted = lift_features(&grouped, &graph)?;

  let mut writer = CsvAnnotationWriter::new(&output, b',')?;
  writer.write_node_annotations(&lifted)?;

  Ok(())
}
