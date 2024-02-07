use crate::commands::build::build_args::{DistanceBackend, PangraphBuildArgs};
use crate::distance::mash::mash_distance::mash_distance;
use crate::distance::mash::minimizer::MinimizersParams;
use crate::io::fasta::{read_many_fasta, FastaRecord};
use crate::io::json::json_write;
use crate::pangraph::pangraph::Pangraph;
use crate::tree::balance::balance;
use crate::tree::neighbor_joining::build_tree_using_neighbor_joining;
use crate::utils::random::get_random_number_generator;
use eyre::Report;
use itertools::Itertools;

pub fn build_run(args: &PangraphBuildArgs) -> Result<(), Report> {
  let PangraphBuildArgs {
    input_fastas,
    len,
    alpha,
    beta,
    circular,
    upper_case,
    sensitivity,
    max_self_map,
    distance_backend,
    alignment_kernel,
    kmer_length,
    seed,
  } = &args;

  let rng = get_random_number_generator(seed);

  let fastas = read_many_fasta(input_fastas)?;

  // TODO: adjust fasta letter case if `upper_case` is set

  // TODO: check for duplicate fasta names

  let pangraph_json = build(fastas, args)?;

  json_write("-", &pangraph_json)?;

  Ok(())
}

pub fn build(fastas: Vec<FastaRecord>, args: &PangraphBuildArgs) -> Result<Vec<Pangraph>, Report> {
  // TODO: check that names are extracted in correct order. Alternatively, consider grouping them with
  //  their respective data entries into a struct, rather than storing names in a separate array.
  let names = fastas.iter().map(|fasta| fasta.seq_name.clone()).collect_vec();

  // Build singleton graphs from input sequences
  let graphs = fastas
    .into_iter()
    .map(|fasta| Pangraph::singleton(fasta, args.circular))
    .collect_vec();

  // Calculate pairwise distances between future guide tree nodes
  let distance = match args.distance_backend {
    DistanceBackend::Native => mash_distance(&graphs, &MinimizersParams::default()),
    DistanceBackend::Mash => {
      // FIXME: what's the difference between Native and Mash?
      unimplemented!("DistanceBackend::Mash");
    }
  };

  // Build guide tree
  let tree = build_tree_using_neighbor_joining(&distance, &names)?;

  // Balance guide tree (to increase available parallelism during parallel traversal?)
  let tree = balance(&tree);

  // TODO: main loop: traverse the tree starting from leaf nodes and and build the graph for the root node

  Ok(graphs)
}
