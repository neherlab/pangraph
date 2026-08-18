use crate::align::alignment_args::check_alignment_backend_available;
use crate::commands::merge::merge_args::PangraphMergeArgs;
use crate::io::json::{JsonPretty, json_write_file};
use crate::make_error;
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_path::PangraphPath;
use crate::pangraph::reconstruct::{check_genome_names, verify_graph_against_graphs};
use color_eyre::owo_colors::{AnsiColors, OwoColorize};
use color_eyre::{Help, SectionExt};
use eyre::{Report, WrapErr};
use log::{info, warn};
use std::collections::BTreeSet;
use std::path::Path;

pub fn merge_run(args: &PangraphMergeArgs) -> Result<(), Report> {
  let left = read_graph(&args.left_graph).wrap_err("When reading the first input graph")?;
  let right = read_graph(&args.right_graph).wrap_err("When reading the second input graph")?;

  merge_cmd_preliminary_checks(args, &left, &right).wrap_err("When performing preliminary checks before merging")?;

  // Block and node ids are derived from genome names, which the checks above established are
  // distinct across the two graphs, so they cannot collide. Path ids are sequential within each
  // graph, so the appended graph is lifted above the first one.
  let right = right
    .renumber_paths(left.path_id_upper_bound())
    .wrap_err("When renumbering the path ids of the second input graph")?;

  // Cheap, and the alternative is `graph_join` panicking on the conflicting key.
  if !right.is_id_disjoint_from(&left) {
    return make_error!(
      "The two input graphs share block or node identifiers, so they cannot be joined. Identifiers are derived from genome names since version 1.4; graphs written by earlier versions derive them from the order of the input sequences instead, and two such graphs collide. Rebuild the input graphs with the current version of pangraph."
    );
  }

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

  if args.verify {
    #[cfg(debug_assertions)]
    merged.sanity_check().wrap_err("When checking the merged graph")?;

    // Compared against the inputs one genome at a time: reconstructing both graphs up front would
    // hold their entire sequence content in memory for the duration of the check.
    verify_graph_against_graphs(&merged, &[&left, &right])
      .wrap_err("When verifying the sequences of the merged graph")?;
    info!(
      "Merged graph reconstructs all {} input genomes exactly",
      merged.paths.len()
    );
  }

  json_write_file(&args.output_json, &merged, JsonPretty(true))?;

  Ok(())
}

/// Reads a pangraph from a JSON file.
///
/// `from_path` already rejects a graph whose ids do not resolve or whose offsets are out of range,
/// in release builds too. The extra `sanity_check` here adds the semantic invariants on top — that
/// node positions tile the genome — which indicate a bug rather than a bad file, and so are only
/// worth paying for in debug builds.
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
  // survives a merge unchanged. Checked unconditionally rather than only under `--verify`: a merge
  // whose genomes cannot be told apart afterwards is not useful either way. Checked here rather
  // than at verification time so that it fails before the expensive merge. `build` enforces the
  // same invariant on its input FASTA records.
  check_genome_names(&[left, right])
    .wrap_err("When checking the genome names of the input graphs")
    .with_section(|| {
      format!("{}\n{}", args.left_graph.display(), args.right_graph.display()).header("Input graphs:")
    })
    .with_section(|| {
      "Genome names must be unique across the two graphs: merging a graph with itself, or re-adding a genome that is already present, is not supported."
        .color(AnsiColors::Cyan)
        .header("Suggestion:")
    })?;

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
