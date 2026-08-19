use crate::align::alignment_args::check_alignment_backend_available;
use crate::commands::build::build_args::PangraphBuildArgs;
use crate::io::fasta::{FastaReader, FastaRecord};
use crate::io::json::{JsonPretty, json_write_file};
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::reconstruct::{GenomeCoverage, check_sequence_names, sequences_by_name, verify_graph_sequences};
use crate::pangraph::strand::Strand::Forward;
use crate::representation::seq::Seq;
use crate::tree::clade::postorder;
use crate::tree::neighbor_joining::build_tree_using_neighbor_joining;
use crate::tree::newick::build_tree_from_newick;
use crate::utils::progress_bar::ProgressBar;
use crate::{make_internal_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;
use std::collections::BTreeMap;

/// Checks that the graph is internally consistent and reconstructs the genomes it should.
fn graph_sanity_checks(
  graph: &Pangraph,
  expected: &BTreeMap<String, Seq>,
  coverage: GenomeCoverage,
) -> Result<(), Report> {
  // check that graph internal structure (blocks, paths, nodes, edits...) is valid
  #[cfg(debug_assertions)]
  graph
    .sanity_check()
    .wrap_err("When performing sanity check on the pangraph")?;

  // Reconstruct the genomes from the graph and compare them with the input sequences. Genomes are
  // matched by name: path ids do not survive a merge, so they cannot pair sequences up.
  verify_graph_sequences(graph, expected, coverage)
    .wrap_err("When comparing reconstructed sequences with original FASTA records")?;

  Ok(())
}

pub fn build_run(args: &PangraphBuildArgs) -> Result<(), Report> {
  let input_fastas = &args.input_fastas;

  let fastas = FastaReader::from_paths(input_fastas)?.read_many()?;

  // TODO: adjust fasta letter case if `upper_case` is set

  check_alignment_backend_available(&args.merge_params)
    .wrap_err("When performing preliminary checks before building the pangraph.")?;

  let pangraph = build(fastas, args, args.verify)?;

  json_write_file(&args.output_json, &pangraph, JsonPretty(true))?;

  Ok(())
}

pub fn build(fastas: Vec<FastaRecord>, args: &PangraphBuildArgs, verify: bool) -> Result<Pangraph, Report> {
  check_sequence_names(&fastas).wrap_err("When checking the names of the input sequences")?;

  // If verification is requested, keep the input sequences, keyed by genome name, to compare them
  // with the sequences reconstructed from the graph. Names were just checked to be unique, so no
  // record is lost to a key collision here.
  let expected = verify.then(|| sequences_by_name(&fastas));

  // Build singleton graphs from input sequences
  // TODO: initial graphs can potentially be constructed when initializing tree clades. This could avoid a lot of boilerplate code.
  let n_paths = fastas.len();
  let graphs = fastas
    .into_iter()
    .map(|fasta| Pangraph::singleton(fasta, Forward, args.circular)) // FIXME: strand hardcoded
    .collect_vec();

  // Build guide tree, or load it from a user-supplied Newick file.
  let tree = match &args.guide_tree {
    Some(path) => build_tree_from_newick(path, graphs)
      .wrap_err_with(|| format!("When loading guide tree from '{}'", path.display()))?,
    None => build_tree_using_neighbor_joining(graphs)?,
  };

  // Log the guide tree topology in Newick format (visible at `info` log level and above).
  info!("Guide tree (newick): {}", tree.read().to_newick());

  // Instantiate the progress bar
  let pb = ProgressBar::new(n_paths - 1, args.no_progress_bar)?;

  // Main loop: traverse the tree starting from leaf nodes and build the graphs bottom-up all the way to the root node.
  // The graph of the root node is the graph we are looking for.
  postorder(&tree, |clade| {
    match (&clade.left, &clade.right) {
      (None, None) => {
        // Case: leaf node. Action: nothing to do.
        Ok(())
      },
      (Some(left), Some(right)) => {
        // Case: internal node with two children. Action: produce graph for this node based on the graphs of its children.
        // Assumption: Child nodes are assumed to be already visited at this point.
        match (&left.read().data, &right.read().data) {
          (Some(left), Some(right)) => {
            info!(
              "=== Graph merging start:     clades sizes {} + {}",
              left.paths.len(),
              right.paths.len()
            );

            clade.data = Some(merge_graphs(left, right, &args.merge_params).wrap_err("When merging graphs")?);

            // increase progress bar
            pb.inc(1);

            info!(
              "=== Graph merging completed: clades sizes {} + {} -> {}",
              left.paths.len(),
              right.paths.len(),
              clade.data.as_ref().unwrap().paths.len()
            );

            // perform checks only in debug mode and if requested. An intermediate graph holds
            // only the genomes of its own clade, hence `Partial`.
            #[cfg(debug_assertions)]
            if let Some(expected) = expected.as_ref() {
              graph_sanity_checks(clade.data.as_ref().unwrap(), expected, GenomeCoverage::Partial)
                .wrap_err("When performing sanity checks on the merged graph")?;
            }

            Ok(())
          },
          _ => {
            make_internal_error!("Found internal clade with two children, of which one or both have no graph attached.")
          },
        }
      },
      (None, Some(_)) | (Some(_), None) => {
        // Case: internal node with one child. Action: ???
        unimplemented!("What to do if there's only one child?");
      },
    }
  })
  .wrap_err("When traversing guide tree")?;

  // Finish progress bar
  pb.finish_with_message("Graph merging completed");

  let graph = tree
    .write()
    .data
    .take()
    .ok_or_else(|| make_internal_report!("Root clade of the guide tree contains no graph after graph alignment"))?;

  // verify the final graph if requested. It must hold every input genome, hence `Complete`.
  if let Some(expected) = &expected {
    graph_sanity_checks(&graph, expected, GenomeCoverage::Complete)
      .wrap_err("When performing sanity checks on the final pangraph")?;
    info!("Pangraph reconstructs all {} input genomes exactly", expected.len());
  }

  Ok(graph)
}
