use crate::commands::generate::generate_args::PangraphGenerateArgs;
use crate::io::fasta::{read_many_fasta, FastaWriter};
use crate::utils::random::get_random_number_generator;
use eyre::Report;

pub fn generate_run(args: &PangraphGenerateArgs) -> Result<(), Report> {
  let PangraphGenerateArgs {
    input_fastas,
    snp_rate,
    hgt_rate,
    sigma_pre,
    delete_rate,
    invert_rate,
    graph_output,
    time,
    seed,
  } = &args;

  let rng = get_random_number_generator(seed);

  let fastas = read_many_fasta(input_fastas)?;

  let mut fasta_writer = FastaWriter::from_path("-")?;
  for fasta in fastas {
    fasta_writer.write(&fasta.seq_name, &fasta.seq)?;
  }

  Ok(())
}
