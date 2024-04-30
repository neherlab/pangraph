use clap::{Parser, ValueHint};
use ctor::ctor;
use eyre::Report;
use itertools::Itertools;
use pangraph::align::alignment_args::AlignmentArgs;
use pangraph::align::mmseqs::align_with_mmseqs::align_with_mmseqs;
use pangraph::io::fasta::read_many_fasta;
use pangraph::utils::global_init::global_init;
use std::path::PathBuf;

#[ctor]
fn init() {
  global_init();
}

#[derive(Parser, Debug)]
struct Args {
  #[clap(value_hint = ValueHint::FilePath)]
  pub input_query_fastas: Vec<PathBuf>,

  #[clap(flatten)]
  pub params: AlignmentArgs,
}

fn main() -> Result<(), Report> {
  let Args {
    input_query_fastas,
    params,
  } = Args::parse();

  let qrys = read_many_fasta(&input_query_fastas)?
    .into_iter()
    .map(|r| r.seq)
    .collect_vec();

  let result = align_with_mmseqs(&qrys, &qrys, &params)?;
  println!("{:#?}", &result);

  Ok(())
}
