use crate::commands::build::build_args::PangraphBuildArgs;
use crate::graph::pangraph::Pangraph;
use crate::io::fasta::{read_many_fasta, FastaRecord};
use crate::io::json::json_write;
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

  let pangraph_json = build(&fastas, args)?;

  json_write("-", &pangraph_json)?;

  Ok(())
}

pub fn build(fastas: &[FastaRecord], args: &PangraphBuildArgs) -> Result<Vec<Pangraph>, Report> {
  // Build singleton graphs from input sequences
  let graphs = fastas
    .iter()
    .map(|fasta| Pangraph::singleton(fasta, args.circular))
    .collect_vec();

  // TODO: build guide tree

  // TODO: balance guide tree

  // TODO: main loop: traverse the tree starting from leaf nodes and and build the graph for the root node

  Ok(graphs)
}
