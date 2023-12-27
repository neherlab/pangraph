use ctor::ctor;
use eyre::Report;
use log::info;
use pangraph::commands::build::build_run::build_run;
use pangraph::commands::export::export_run::export_run;
use pangraph::commands::generate::generate_run::generate_run;
use pangraph::commands::marginalize::marginalize_run::marginalize_run;
use pangraph::commands::polish::polish_run::polish_run;
use pangraph::commands::root_args::{generate_shell_completions, parse_cli_args, PangraphCommands};
use pangraph::utils::global_init::global_init;

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

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
    PangraphCommands::Generate(args) => generate_run(&args),
    PangraphCommands::Marginalize(args) => marginalize_run(&args),
    PangraphCommands::Polish(args) => polish_run(&args),
    PangraphCommands::Completions { shell } => generate_shell_completions(&shell),
  }
}
