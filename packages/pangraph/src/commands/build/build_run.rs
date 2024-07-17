use crate::commands::build::build_args::PangraphBuildArgs;
use crate::commands::reconstruct::reconstruct_run::{compare_sequences, reconstruct};
use crate::graph::traversal::GraphTraversalContinuation;
use crate::io::fasta::{read_many_fasta, FastaRecord};
use crate::io::json::{json_stringify, json_write};
use crate::make_internal_report;
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::strand::Strand::Forward;
use crate::tree::neighbor_joining::build_tree_using_neighbor_joining;
use crate::utils::random::get_random_number_generator;
use eyre::{Report, WrapErr};
use itertools::Itertools;

pub fn build_run(args: &PangraphBuildArgs) -> Result<(), Report> {
  let PangraphBuildArgs { input_fastas, seed, .. } = &args;

  let rng = get_random_number_generator(seed);

  let fastas = read_many_fasta(input_fastas)?;

  // TODO: adjust fasta letter case if `upper_case` is set
  // TODO: check for duplicate fasta names

  let pangraph = if args.verify {
    let pangraph = build(fastas.clone(), args)?;
    {
      let mut results = reconstruct(&pangraph);
      results.try_for_each(|actual| -> Result<(), Report> {
        let actual = actual?;
        let expected = &fastas[actual.index];
        compare_sequences(expected, &actual);
        Ok(())
      })?;
    }
    Ok(pangraph)
  } else {
    build(fastas, args)
  }?;

  json_write(&args.output_json, &pangraph)?;

  Ok(())
}

pub fn build(fastas: Vec<FastaRecord>, args: &PangraphBuildArgs) -> Result<Pangraph, Report> {
  // Build singleton graphs from input sequences
  // TODO: initial graphs can potentially be constructed when initializing tree clades. This could avoid a lot of boilerplate code.
  let graphs = fastas
    .into_iter()
    .map(|fasta| Pangraph::singleton(fasta, Forward, args.circular)) // FIXME: strand hardcoded
    .collect_vec();

  // Build guide tree
  let tree = build_tree_using_neighbor_joining(graphs, args)?;

  // Main loop: traverse the tree starting from leaf nodes and build the graphs bottom-up all the way to the root node.
  // The graph of the root node is the graph we are looking for.
  tree.traverse_backward(|node| {
    let mut clade = node.payload;

    if node.children.is_empty() {
      return GraphTraversalContinuation::Continue;
    }

    let left = &node.children[0].0.read_arc().pangraph;
    let right = &node.children[1].0.read_arc().pangraph;

    match (left, right) {
      (None, None) => {
        // Case: leaf node. Action: nothing to do.
      }
      (Some(left), Some(right)) => {
        // Case: internal node with two children. Action: produce graph for this node based on the graphs of its children.
        // Assumption: Child nodes are assumed to be already visited at this point.
        clade.pangraph = Some(merge_graphs(left, right, args).wrap_err("When merging graphs").unwrap());
      }
      (None, Some(child)) | (Some(child), None) => {
        // Case: internal node with one child. Action: ???
        unimplemented!("What to do if there's only one child?");
      }
    }
    GraphTraversalContinuation::Continue
  });

  let graph = tree
    .get_exactly_one_root()?
    .write()
    .payload()
    .write()
    .pangraph
    .take()
    .ok_or_else(|| make_internal_report!("Root clade of the guide tree contains no graph after graph alignment"))?;

  Ok(graph)
}
