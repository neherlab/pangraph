use crate::commands::simplify::simplify_args::PangraphSimplifyArgs;
use crate::pangraph::pangraph::Pangraph;
use crate::utils::random::get_random_number_generator;
use eyre::Report;

pub fn simplify_run(args: &PangraphSimplifyArgs) -> Result<(), Report> {
  let PangraphSimplifyArgs {
    input_aln,
    output_path,
    strains,
    seed,
  } = &args;

  let rng = get_random_number_generator(seed);

  let msa_json = Pangraph::from_path(input_aln)?;

  Ok(())
}
