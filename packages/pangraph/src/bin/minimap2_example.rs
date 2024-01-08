use ctor::ctor;
use eyre::Report;
use pangraph::align::minimap2::align_with_minimap2::align_with_minimap2;
use pangraph::utils::global_init::global_init;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let results = align_with_minimap2("data/smol_01_ref.fa", "data/smol_01.fa")?;
  dbg!(&results);
  Ok(())
}
