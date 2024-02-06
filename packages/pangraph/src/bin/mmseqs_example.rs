use ctor::ctor;
use eyre::Report;
use pangraph::align::mmseqs::align_with_mmseqs::{align_with_mmseqs, MmseqsParams};
use pangraph::io::fasta::read_one_fasta;
use pangraph::utils::global_init::global_init;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let ref_seq = read_one_fasta("data/smol_01_ref.fa")?;
  let qry_seq = read_one_fasta("data/smol_01.fa")?;

  let results = align_with_mmseqs(&ref_seq, &qry_seq, &MmseqsParams::default())?;
  dbg!(&results);
  Ok(())
}
