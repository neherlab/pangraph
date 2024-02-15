use clap::{ArgEnum, Parser, ValueHint};
use ctor::ctor;
use eyre::Report;
use pangraph::align::mmseqs::align_with_mmseqs::{align_with_mmseqs, MmseqsParams};
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

  #[clap(flatten)]
  pub params: MmseqsParams,
}

fn main() -> Result<(), Report> {
  let Args {
    input_ref_fasta: input_ref_fastas,
    input_query_fastas,
    params,
  } = Args::parse();

  let ref_record = read_one_fasta(input_ref_fastas)?;

  let qry_records = read_many_fasta(&input_query_fastas)?;

  for qry_record in qry_records {
    let result = align_with_mmseqs(&ref_record.seq, &qry_record.seq, &params)?;
    dbg!(&result);
  }

  Ok(())
}
