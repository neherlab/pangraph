use crate::commands::build::build_run::build_run;
use crate::commands::export::export_run::export_run;
use crate::commands::marginalize::marginalize_run::marginalize_run;
use crate::commands::md_help::print_help_markdown::print_help_markdown;
use crate::commands::reconstruct::reconstruct_run::reconstruct_run;
use crate::commands::root_args::{generate_shell_completions, parse_cli_args, PangraphCommands};
use crate::commands::schema::generate_schema::generate_schema;
use eyre::Report;
use log::info;

pub fn pangraph_main() -> Result<(), Report> {
  let args = parse_cli_args()?;

  info!("{:#?}", &args);

  rayon::ThreadPoolBuilder::new().num_threads(args.jobs).build_global()?;

  match args.command {
    PangraphCommands::Build(args) => build_run(&args),
    PangraphCommands::Export(args) => export_run(&args),
    PangraphCommands::Marginalize(args) => marginalize_run(&args),
    PangraphCommands::Reconstruct(args) => reconstruct_run(&args),
    PangraphCommands::Schema(args) => generate_schema(&args),
    PangraphCommands::HelpMarkdown => print_help_markdown(),
    PangraphCommands::Completions { shell } => generate_shell_completions(&shell),
  }
}
