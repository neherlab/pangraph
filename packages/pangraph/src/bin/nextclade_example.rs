use clap::{Parser, ValueHint};
use ctor::ctor;
use eyre::Report;
use pangraph::align::map_variations::BandParameters;
use pangraph::align::nextclade::align_with_nextclade::{NextalignParams, align_with_nextclade};
use pangraph::io::fasta::FastaReader;
use pangraph::utils::global_init::{global_init, setup_logger};
use std::path::PathBuf;

#[ctor]
fn init() {
  global_init();
  setup_logger(log::LevelFilter::Warn);
}

#[derive(Parser, Debug)]
struct Args {
  #[clap(long, short = 'r')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub input_ref_fasta: PathBuf,

  #[clap(value_hint = ValueHint::FilePath)]
  pub input_query_fastas: Vec<PathBuf>,

  #[clap(long, short = 's')]
  pub mean_shift: i32,

  #[clap(long, short = 'b')]
  pub bandwidth: usize,
}

fn main() -> Result<(), Report> {
  let Args {
    input_ref_fasta,
    input_query_fastas,
    mean_shift,
    bandwidth,
  } = Args::parse();

  let ref_record = FastaReader::from_path(input_ref_fasta)?.read_one()?;
  let params = NextalignParams {
    max_alignment_attempts: 1,
    ..Default::default()
  };
  let band_params = BandParameters::new(mean_shift, bandwidth);

  let qry_records = FastaReader::from_paths(&input_query_fastas)?.read_many()?;
  for qry_record in qry_records {
    let result = align_with_nextclade(&ref_record.seq, &qry_record.seq, band_params, &params)?;
    println!("{:#?}", &result);
  }

  Ok(())
}
