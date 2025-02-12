use ctor::ctor;
use eyre::Report;
use pangraph::commands::main::pangraph_main;
use pangraph::utils::global_init::global_init;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  pangraph_main()
}
