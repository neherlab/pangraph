use clap::{ArgEnum, Parser, ValueHint};
use ctor::ctor;
use eyre::Report;
use pangraph::align::nextclade::align_with_nextclade::{align_with_nextclade, NextalignParams};
use pangraph::io::fasta::{read_many_fasta, read_one_fasta};
use pangraph::o;
use pangraph::utils::global_init::global_init;
use std::path::PathBuf;

#[ctor]
fn init() {
  global_init();
}

#[derive(Parser, Debug)]
struct Args {
  #[clap(long, short = 'r')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub input_ref_fasta: PathBuf,

  #[clap(value_hint = ValueHint::FilePath)]
  pub input_query_fastas: Vec<PathBuf>,
}

fn main() -> Result<(), Report> {
  let Args {
    input_ref_fasta,
    input_query_fastas,
  } = Args::parse();

  let ref_record = read_one_fasta(input_ref_fasta)?;

  let qry_records = read_many_fasta(&input_query_fastas)?;

  for qry_record in qry_records {
    let result = align_with_nextclade(&ref_record.seq, &qry_record.seq, &NextalignParams::default())?;
    dbg!(&result);
  }

  Ok(())
}
