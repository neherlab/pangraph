use crate::align::alignment_args::check_alignment_backend_available;
use crate::commands::merge::merge_args::PangraphMergeArgs;
use crate::io::json::{JsonPretty, json_write_file};
use crate::make_error;
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_path::PangraphPath;
use crate::pangraph::reconstruct::reconstruct;
use crate::representation::seq::Seq;
use crate::utils::collections::find_duplicates;
use eyre::{Report, WrapErr};
use itertools::Itertools;
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
    .then(|| expected_sequences(&left, &right))
    .transpose()
    .wrap_err("When reconstructing the sequences of the input graphs")?;

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
    verify_merged_sequences(&merged, &expected).wrap_err("When verifying the sequences of the merged graph")?;
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

  if args.verify {
    for (graph, filepath) in [(left, &args.left_graph), (right, &args.right_graph)] {
      if graph.path_names().any(|name| name.is_none()) {
        return make_error!(
          "Graph '{}' contains genomes without a name, which cannot be verified: verification matches genomes by name. Re-run without `--verify`.",
          filepath.display()
        );
      }
    }
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
fn expected_sequences(left: &Pangraph, right: &Pangraph) -> Result<BTreeMap<String, Seq>, Report> {
  let mut expected = BTreeMap::new();
  for graph in [left, right] {
    for record in reconstruct(graph) {
      let record = record?;
      expected.insert(record.seq_name, record.seq);
    }
  }
  Ok(expected)
}

/// Checks that the merged graph reconstructs exactly the genomes of the input graphs.
/// Genomes are matched by name: path ids are renumbered by the merger, and the order in which
/// genomes are reconstructed is therefore not the order of either input graph.
fn verify_merged_sequences(merged: &Pangraph, expected: &BTreeMap<String, Seq>) -> Result<(), Report> {
  #[cfg(debug_assertions)]
  merged.sanity_check().wrap_err("When checking the merged graph")?;

  let mut remaining: BTreeSet<&String> = expected.keys().collect();

  for record in reconstruct(merged) {
    let record = record?;
    let Some(expected_seq) = expected.get(&record.seq_name) else {
      return make_error!(
        "Merged graph contains genome '{}', which is not present in either input graph",
        record.seq_name
      );
    };

    if record.seq != *expected_seq {
      return make_error!(
        "Sequence mismatch for genome '{}': expected length {} but got {}",
        record.seq_name,
        expected_seq.len(),
        record.seq.len()
      );
    }

    remaining.remove(&record.seq_name);
  }

  if !remaining.is_empty() {
    return make_error!(
      "Merged graph is missing genomes from the input graphs: [{}]",
      remaining.into_iter().sorted().join(", ")
    );
  }

  Ok(())
}
