#![allow(unused_qualifications)]

use crate::commands::build::build_args::PangraphBuildArgs;
use crate::commands::export::export_args::PangraphExportArgs;
use crate::commands::reconstruct::reconstruct_args::PangraphReconstructArgs;
use crate::commands::schema::generate_schema::PangraphGenerateSchemaArgs;
use crate::commands::simplify::simplify_args::PangraphSimplifyArgs;
use crate::commands::verbosity::Verbosity;
use crate::utils::global_init::setup_logger;
use clap::builder::styling;
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::{Shell, generate};
use clap_complete_fig::Fig;
use eyre::{Report, eyre};
use num_cpus;
use std::fmt::Debug;
use std::io;

const SHELLS: &[&str] = &["bash", "elvish", "fish", "fig", "powershell", "zsh"];

fn styles() -> styling::Styles {
  styling::Styles::styled()
    .header(styling::AnsiColor::Green.on_default() | styling::Effects::BOLD)
    .usage(styling::AnsiColor::Green.on_default() | styling::Effects::BOLD)
    .literal(styling::AnsiColor::Blue.on_default() | styling::Effects::BOLD)
    .placeholder(styling::AnsiColor::Cyan.on_default())
}

#[derive(Parser, Debug)]
#[clap(name = "pangraph")]
#[clap(author, version)]
#[clap(verbatim_doc_comment)]
#[clap(styles = styles())]
/// Bioinformatic toolkit to align large sets of closely related genomes into a graph data structure.
///
/// Finds homology amongst large collections of closely related genomes. The core of the algorithm partitions each genome into pancontigs (also called blocks) that represent a sequence interval related by vertical descent. Each genome is then an ordered walk along pancontigs. The collection of all genomes form a graph that captures all observed structural diversity. The tool is useful to study structural variations in the genome, perform comparative studies of genome gain, loss, and rearrangement dynamics; or simply to compress many related genomes.
///
///
/// Publication: "PanGraph: scalable bacterial pan-genome graph construction." Nicholas Noll, Marco Molari, Richard Neher. Microbial Genomics 9.6 (2023): 001034.; doi: https://doi.org/10.1099/mgen.0.001034
///
/// Documentation: https://pangraph.readthedocs.io/en/stable/
///
/// Source code: https://github.com/neherlab/pangraph
///
/// Questions, ideas, bug reports: https://github.com/neherlab/pangraph/issues
pub struct PangraphArgs {
  #[clap(subcommand)]
  pub command: PangraphCommands,

  /// Make output more quiet or more verbose
  #[clap(flatten, next_help_heading = "Verbosity")]
  pub verbosity: Verbosity,

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
  Export {
    #[clap(subcommand)]
    args: PangraphExportArgs,
  },

  /// Generates a simplified graph that only contains a subset of the input genomes.
  Simplify(PangraphSimplifyArgs),

  /// Reconstruct all input fasta sequences from graph
  Reconstruct(PangraphReconstructArgs),

  /// Generate JSON schema for Pangraph file format
  Schema(PangraphGenerateSchemaArgs),

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
    #[clap(value_name = "SHELL", default_value_t = String::from("bash"), value_parser = SHELLS.to_vec())]
    shell: String,
  },

  /// Print command-line reference documentation in Markdown format
  HelpMarkdown,
}

pub fn generate_shell_completions(shell: &str) -> Result<(), Report> {
  let mut command = PangraphArgs::command();

  if shell.to_lowercase() == "fig" {
    generate(Fig, &mut command, "pangraph", &mut io::stdout());
    return Ok(());
  }

  let generator = <Shell as ValueEnum>::from_str(&shell.to_lowercase(), true)
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
