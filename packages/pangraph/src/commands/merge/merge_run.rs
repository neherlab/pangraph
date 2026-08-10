use crate::align::alignment_args::check_alignment_backend_available;
use crate::commands::merge::merge_args::PangraphMergeArgs;
use crate::io::json::{JsonPretty, json_write_file};
use crate::make_error;
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_path::PangraphPath;
use crate::pangraph::reconstruct::{GenomeCoverage, reconstruct_by_name, verify_graph_sequences};
use crate::representation::seq::Seq;
use crate::utils::collections::find_duplicates;
use color_eyre::owo_colors::{AnsiColors, OwoColorize};
use color_eyre::{Help, SectionExt};
use eyre::{Report, WrapErr};
use log::{info, warn};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub fn merge_run(args: &PangraphMergeArgs) -> Result<(), Report> {
  let left = read_graph(&args.left_graph).wrap_err("When reading the first input graph")?;
  let right = read_graph(&args.right_graph).wrap_err("When reading the second input graph")?;

  merge_cmd_preliminary_checks(args, &left, &right).wrap_err("When performing preliminary checks before merging")?;

  // The two graphs were built independently, so their identifiers almost certainly collide.
  // Namespace the second graph before joining them.
  let right = right
    .make_disjoint_from(&left)
    .wrap_err("When making the identifiers of the two input graphs disjoint")?;

  // Reconstruct the expected genomes from the (relabeled) inputs, before they are consumed by the
  // merger. Keyed by path name: neither path ids nor record order survive a merge.
  let expected = args
    .verify
    .then(|| expected_sequences(args, &left, &right))
    .transpose()?;

  info!(
    "=== Graph merging start:     graph sizes {} + {}",
    left.paths.len(),
    right.paths.len()
  );

  let merged = merge_graphs(&left, &right, &args.merge_params).wrap_err("When merging graphs")?;

  info!(
    "=== Graph merging completed: graph sizes {} + {} -> {} paths, {} blocks",
    left.paths.len(),
    right.paths.len(),
    merged.paths.len(),
    merged.blocks.len()
  );

  if let Some(expected) = expected {
    #[cfg(debug_assertions)]
    merged.sanity_check().wrap_err("When checking the merged graph")?;

    verify_graph_sequences(&merged, &expected, GenomeCoverage::Complete)
      .wrap_err("When verifying the sequences of the merged graph")?;
    info!("Merged graph reconstructs all {} input genomes exactly", expected.len());
  }

  json_write_file(&args.output_json, &merged, JsonPretty(true))?;

  Ok(())
}

/// Reads a pangraph from a JSON file, and checks its internal consistency in debug builds.
fn read_graph(filepath: &Path) -> Result<Pangraph, Report> {
  let graph = Pangraph::from_path(&Some(filepath))?;

  #[cfg(debug_assertions)]
  graph
    .sanity_check()
    .wrap_err_with(|| format!("When performing sanity check on graph '{}'", filepath.display()))?;

  Ok(graph)
}

fn merge_cmd_preliminary_checks(args: &PangraphMergeArgs, left: &Pangraph, right: &Pangraph) -> Result<(), Report> {
  check_alignment_backend_available(&args.merge_params)?;

  for (graph, filepath) in [(left, &args.left_graph), (right, &args.right_graph)] {
    if graph.paths.is_empty() {
      return make_error!("Input graph '{}' contains no genomes", filepath.display());
    }
  }

  // Genome names are the identity of a genome throughout pangraph, and the only handle that
  // survives a merge unchanged. `build` enforces the same invariant on its input FASTA records.
  let duplicates = find_duplicates([left, right].into_iter().flat_map(|graph| graph.path_names().flatten()));
  if !duplicates.is_empty() {
    return make_error!(
      "Duplicate genome names found in the input graphs: [{}]. Genome names must be unique across the two graphs: merging a graph with itself, or re-adding a genome that is already present, is not supported.",
      duplicates.join(", ")
    );
  }

  // Circularity is a per-path property, so mixing is structurally fine. It is however most often a
  // mistake, since `build --circular` applies to all genomes of a graph at once.
  if circularity(left) != circularity(right) {
    warn!(
      "The two input graphs disagree on whether their genomes are circular. This is allowed, but check that both graphs were built with a consistent `--circular` setting."
    );
  }

  Ok(())
}

/// Returns the set of circularity flags used by the paths of a graph.
fn circularity(graph: &Pangraph) -> BTreeSet<bool> {
  graph.paths().map(PangraphPath::circular).collect()
}

/// Reconstructs the genomes of both input graphs, keyed by genome name.
///
/// Cross-graph name collisions are already rejected by `merge_cmd_preliminary_checks`, so the two
/// sets cannot overwrite each other here.
fn expected_sequences(
  args: &PangraphMergeArgs,
  left: &Pangraph,
  right: &Pangraph,
) -> Result<BTreeMap<String, Seq>, Report> {
  let mut expected = BTreeMap::new();
  for (graph, filepath) in [(left, &args.left_graph), (right, &args.right_graph)] {
    let genomes = reconstruct_by_name(graph)
      .wrap_err_with(|| format!("When reconstructing the genomes of graph '{}'", filepath.display()))
      .with_section(|| {
        "Verification matches genomes by name. Re-run without `--verify` to skip it."
          .color(AnsiColors::Cyan)
          .header("Suggestion:")
      })?;
    expected.extend(genomes);
  }
  Ok(expected)
}
