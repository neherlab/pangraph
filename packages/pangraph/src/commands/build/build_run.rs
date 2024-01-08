use crate::commands::build::build_args::PangraphBuildArgs;
use crate::graph::pangraph::Pangraph;
use crate::io::fasta::{read_many_fasta, FastaRecord};
use crate::io::json::json_write;
use crate::utils::random::get_random_number_generator;
use eyre::Report;

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
  let pangraph_json = build(&fastas, args)?;
  json_write("-", &pangraph_json)?;

  Ok(())
}

pub fn build(fastas: &[FastaRecord], args: &PangraphBuildArgs) -> Result<Pangraph, Report> {
  // TODO: create proper pangraph JSON ;)
  let pangraph_json = Pangraph {
    paths: vec![],
    blocks: vec![],
  };
  Ok(pangraph_json)
}
