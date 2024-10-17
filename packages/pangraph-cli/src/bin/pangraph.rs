use ctor::ctor;
use eyre::Report;
use log::info;
use pangraph::commands::build::build_run::build_run;
use pangraph::commands::export::export_run::export_run;
use pangraph::commands::marginalize::marginalize_run::marginalize_run;
use pangraph::commands::reconstruct::reconstruct_run::reconstruct_run;
use pangraph::commands::root_args::{generate_shell_completions, parse_cli_args, PangraphCommands};
use pangraph::commands::schema::generate_schema::generate_schema;
use pangraph::utils::global_init::global_init;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let args = parse_cli_args()?;

  info!("{:#?}", &args);

  rayon::ThreadPoolBuilder::new().num_threads(args.jobs).build_global()?;

  match args.command {
    PangraphCommands::Build(args) => build_run(&args),
    PangraphCommands::Export(args) => export_run(&args),
    PangraphCommands::Marginalize(args) => marginalize_run(&args),
    PangraphCommands::Reconstruct(args) => reconstruct_run(&args),
    PangraphCommands::Schema(args) => generate_schema(&args),
    PangraphCommands::Completions { shell } => generate_shell_completions(&shell),
  }
}
