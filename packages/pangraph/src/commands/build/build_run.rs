use crate::align::align_graphs::align_graphs;
use crate::commands::build::build_args::PangraphBuildArgs;
use crate::io::fasta::{read_many_fasta, FastaRecord};
use crate::io::json::json_write;
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
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
  // Build singleton graphs from input sequences
  // TODO: initial graphs can potentially be constructed when initializing tree clades. This could avoid a lot of boilerplate code.
  let graphs = fastas
    .into_iter()
    .map(|fasta| Pangraph::singleton(fasta, true, args.circular)) // FIXME: strand hardcoded
    .collect_vec();

  // Build guide tree
  let tree = build_tree_using_neighbor_joining(graphs, args)?;

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
        if let (Some(left), Some(right)) = (&left.read().data, &right.read().data) {
          clade.data = Some(merge_graphs(left, right, args)?);
          Ok(())
        } else {
          make_internal_error!("Found internal clade with two children, of which one or both have no graph attached.")
        }
      }
      (None, Some(child)) | (Some(child), None) => {
        // Case: internal node with one child. Action: ???
        unimplemented!("What to do if there's only one child?");
      }
    }
  })
  .into_iter()
  .collect::<Result<Vec<_>, Report>>()
  .wrap_err("When traversing guide tree")?;

  let graph = tree
    .write()
    .data
    .take()
    .ok_or_else(|| make_internal_report!("Root clade of the guide tree contains no graph after graph alignment"))?;

  Ok(graph)
}
