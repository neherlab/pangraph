use crate::align::align_graphs::align_graphs;
use crate::commands::build::build_args::{DistanceBackend, PangraphBuildArgs};
use crate::distance::mash::mash_distance::mash_distance;
use crate::distance::mash::minimizer::MinimizersParams;
use crate::io::fasta::{read_many_fasta, FastaRecord};
use crate::io::json::json_write;
use crate::pangraph::pangraph::Pangraph;
use crate::tree::balance::balance;
use crate::tree::clade::postorder;
use crate::tree::neighbor_joining::build_tree_using_neighbor_joining;
use crate::utils::random::get_random_number_generator;
use crate::{make_internal_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;

pub fn build_run(args: &PangraphBuildArgs) -> Result<(), Report> {
  let PangraphBuildArgs { input_fastas, seed, .. } = &args;

  let rng = get_random_number_generator(seed);

  let fastas = read_many_fasta(input_fastas)?;

  // TODO: adjust fasta letter case if `upper_case` is set

  // TODO: check for duplicate fasta names

  let pangraph_json = build(fastas, args)?;

  json_write("-", &pangraph_json)?;

  Ok(())
}

pub fn build(fastas: Vec<FastaRecord>, args: &PangraphBuildArgs) -> Result<Pangraph, Report> {
  // TODO: check that names are extracted in correct order. Alternatively, consider grouping them with
  //  their respective data entries into a struct, rather than storing names in a separate array.
  let names = fastas.iter().map(|fasta| fasta.seq_name.clone()).collect_vec();

  // Build singleton graphs from input sequences
  // TODO: initial graphs can potentially be constructed when initializing tree clades. This could avoid a lot of boilerplate code.
  let graphs = fastas
    .into_iter()
    .map(|fasta| Pangraph::singleton(fasta, args.circular))
    .collect_vec();

  // Calculate pairwise distances between future guide tree nodes
  let distance = match args.distance_backend {
    // TODO: this function only needs sequences, and not graphs
    DistanceBackend::Native => mash_distance(&graphs, &MinimizersParams::default()),
    DistanceBackend::Mash => {
      // FIXME: what's the difference between Native and Mash?
      unimplemented!("DistanceBackend::Mash");
    }
  };

  // Build guide tree
  let tree = build_tree_using_neighbor_joining(&distance, &names, &graphs)?;

  // Balance guide tree (to increase available parallelism during parallel traversal?)
  let tree = balance(&tree);

  // Main loop: traverse the tree starting from leaf nodes and build the graphs bottom-up all the way to the root node.
  // The graph of the root node is the graph we are looking for.
  postorder(&tree, |clade| {
    match (&clade.left, &clade.right) {
      (None, None) => {
        // Case: leaf node. Action: nothing to do.
        Ok(())
      }
      (Some(left), Some(right)) => {
        // Case: internal node with two children. Action: produce graph for this node based on the graphs of its children.
        // Assumption: Child nodes are assumed to be already visited at this point.
        if let (Some(left), Some(right)) = (&left.read().graph, &right.read().graph) {
          clade.graph = Some(align_graphs(left, right, args)?);
          Ok(())
        } else {
          make_internal_error!("Found internal clade with two children, of which one or both have no graph attached.")
        }
      }
      (None, Some(child)) | (Some(child), None) => {
        // Case: internal node with one child. Action: ???
        todo!("What to do if there's only one child?");
      }
    }
  })
  .into_iter()
  .collect::<Result<Vec<_>, Report>>()
  .wrap_err("When traversing guide tree")?;

  let graph = tree
    .write()
    .graph
    .take()
    .ok_or_else(|| make_internal_report!("Root clade of the guide tree contains no graph after graph alignment"))?;

  Ok(graph)
}
