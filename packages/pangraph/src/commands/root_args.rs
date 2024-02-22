use crate::commands::build::build_args::PangraphBuildArgs;
use crate::commands::export::export_args::PangraphExportArgs;
use crate::commands::generate::generate_args::PangraphGenerateArgs;
use crate::commands::marginalize::marginalize_args::PangraphMarginalizeArgs;
use crate::commands::polish::polish_args::PangraphPolishArgs;
use crate::commands::verbosity::{Verbosity, WarnLevel};
use crate::utils::global_init::setup_logger;
use clap::{AppSettings, ArgEnum, CommandFactory, Parser, Subcommand};
use clap_complete::{generate, Generator, Shell};
use clap_complete_fig::Fig;
use eyre::{eyre, Report};
use itertools::Itertools;
use lazy_static::lazy_static;
use num_cpus;
use std::fmt::Debug;
use std::io;
use std::str::FromStr;

lazy_static! {
  static ref SHELLS: &'static [&'static str] = &["bash", "elvish", "fish", "fig", "powershell", "zsh"];
  static ref VERBOSITIES: &'static [&'static str] = &["off", "error", "warn", "info", "debug", "trace"];
}

#[derive(Parser, Debug)]
#[clap(name = "pangraph", trailing_var_arg = true)]
#[clap(author, version)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(verbatim_doc_comment)]
/// Bioinformatic toolkit to align large sets of closely related genomes into a graph data structure.
///
/// Finds homology amongst large collections of closely related genomes. The core of the algorithm partitions each genome into pancontigs that represent a sequence interval related by vertical descent. Each genome is then an ordered walk along pancontigs; the collection of all genomes form a graph that captures all observed structural diversity. The tool useful to parsimoniously infer horizontal gene transfer events within a community; perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.
///
///
/// Publication: "PanGraph: scalable bacterial pan-genome graph construction. Nicholas Noll, Marco Molari, Richard Neher. bioRxiv 2022.02.24.481757; doi: https://doi.org/10.1101/2022.02.24.481757"
///
/// Documentation: https://pangraph.readthedocs.io/en/stable/
///
/// Source code:https://github.com/neherlab/pangraph
///
/// Questions, ideas, bug reports: https://github.com/neherlab/pangraph/issues
pub struct PangraphArgs {
  #[clap(subcommand)]
  pub command: PangraphCommands,

  /// Make output more quiet or more verbose
  #[clap(flatten, next_help_heading = "  Verbosity")]
  pub verbosity: Verbosity<WarnLevel>,

  /// Number of processing jobs. If not specified, all available CPU threads will be used.
  #[clap(global = true, long, short = 'j', default_value_t = num_cpus::get())]
  pub jobs: usize,
}

#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum PangraphCommands {
  /// Align genomes into a multiple sequence alignment graph
  Build(PangraphBuildArgs),

  /// Export a pangraph to a chosen file format(s)
  Export(PangraphExportArgs),

  /// Output a simulated sequence alignment graph
  Generate(PangraphGenerateArgs),

  /// Compute all pairwise marginalizations of a multiple sequence alignment graph
  Marginalize(PangraphMarginalizeArgs),

  /// Realign pancontigs of multiple sequence alignment graph
  Polish(PangraphPolishArgs),

  /// Generate shell completions.
  ///
  /// This will print the completions file contents to the console. Refer to your shell's documentation on how to install the completions.
  ///
  /// Example for Ubuntu Linux:
  ///
  ///    pangraph completions bash > ~/.local/share/bash-completion/pangraph
  ///
  Completions {
    /// Name of the shell to generate appropriate completions
    #[clap(value_name = "SHELL", default_value_t = String::from("bash"), possible_values(SHELLS.iter()))]
    shell: String,
  },
}

pub fn generate_shell_completions(shell: &str) -> Result<(), Report> {
  let mut command = PangraphArgs::command();

  if shell.to_lowercase() == "fig" {
    generate(Fig, &mut command, "pangraph", &mut io::stdout());
    return Ok(());
  }

  let generator = <Shell as ArgEnum>::from_str(&shell.to_lowercase(), true)
    .map_err(|err| eyre!("{}: Possible values: {}", err, SHELLS.join(", ")))?;

  let bin_name = command.get_name().to_owned();

  generate(generator, &mut command, bin_name, &mut io::stdout());

  Ok(())
}

pub fn parse_cli_args() -> Result<PangraphArgs, Report> {
  let args = PangraphArgs::parse();
  setup_logger(args.verbosity.get_filter_level());
  Ok(args)
}