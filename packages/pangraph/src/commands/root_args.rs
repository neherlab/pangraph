#![allow(unused_qualifications)]

use crate::commands::build::build_args::PangraphBuildArgs;
use crate::commands::export::export_args::PangraphExportArgs;
use crate::commands::merge::merge_args::PangraphMergeArgs;
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
/// Documentation: https://docs.pangraph.org/
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
  // Declared after the `Verbosity` flatten, which opens a help section that would otherwise claim
  // every argument registered after it. A heading set on the argument itself wins over that
  // section, and keeps this option in the default one.
  #[clap(global = true, long, short = 'j', default_value_t = num_cpus::get())]
  #[clap(help_heading = None)]
  pub jobs: usize,
}

#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
pub enum PangraphCommands {
  /// Align genomes into a multiple sequence alignment graph
  Build(PangraphBuildArgs),

  /// Merge two pangenome graphs into a single one
  Merge(PangraphMergeArgs),

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

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::path::PathBuf;

  /// The `Default` impl and the clap `default_value` must agree. They are declared separately, so a
  /// field can easily get one and not the other: an output path then defaults to `""` in code that
  /// builds args programmatically with `..Default::default()`, and writes to a file literally named
  /// `""` instead of to stdout.
  #[rstest]
  fn test_clap_and_rust_defaults_agree_on_output_paths() {
    let stdout = PathBuf::from("-");

    let PangraphCommands::Build(build) = PangraphArgs::parse_from(["pangraph", "build"]).command else {
      panic!("expected the build subcommand");
    };
    assert_eq!(build.output_json, stdout);
    assert_eq!(PangraphBuildArgs::default().output_json, stdout);

    let PangraphCommands::Merge(merge) =
      PangraphArgs::parse_from(["pangraph", "merge", "left.json", "right.json"]).command
    else {
      panic!("expected the merge subcommand");
    };
    assert_eq!(merge.output_json, stdout);
    assert_eq!(PangraphMergeArgs::default().output_json, stdout);

    let PangraphCommands::Reconstruct(reconstruct) = PangraphArgs::parse_from(["pangraph", "reconstruct"]).command
    else {
      panic!("expected the reconstruct subcommand");
    };
    assert_eq!(reconstruct.output_fasta, stdout);
    assert_eq!(PangraphReconstructArgs::default().output_fasta, stdout);
  }

  /// Returns the ids of the arguments of `cmd` filed under the given help section.
  fn args_under(cmd: &clap::Command, heading: &str) -> BTreeSet<String> {
    cmd
      .get_arguments()
      .filter(|arg| arg.get_help_heading() == Some(heading))
      .map(|arg| arg.get_id().to_string())
      .collect()
  }

  fn ids(names: &[&str]) -> BTreeSet<String> {
    names.iter().map(|name| (*name).to_owned()).collect()
  }

  /// A help section is not a property of a group of arguments: clap keeps a single cursor on the
  /// `Command` and stamps it onto each argument as it is registered, and `#[clap(flatten)]` shares
  /// that `Command` with the flattened struct. A flattened group that opens a section therefore
  /// opens it for every argument registered afterwards, including later fields of the *parent*
  /// struct. Declaration order is load-bearing, and getting it wrong is invisible outside `--help`.
  ///
  /// This pins the assignment so that mistake is a test failure. An argument appended after a
  /// flattened group shows up in one of these sets; an argument added before it does not, so there
  /// are no false alarms. The escape hatch, if a trailing argument really is needed, is to set
  /// `#[clap(help_heading = ...)]` on the argument itself — that wins over the cursor.
  #[rstest]
  fn test_arguments_are_filed_under_the_expected_help_section() {
    // `jobs` is declared after the `Verbosity` flatten and used to be swept into it.
    assert_eq!(
      args_under(&PangraphArgs::command(), "Verbosity"),
      ids(&["verbosity", "silent", "verbose", "quiet"])
    );

    // Exactly the options of `GraphMergeParams`, which `build` and `merge` both flatten last.
    let alignment = ids(&[
      "indel_len_threshold",
      "alpha",
      "beta",
      "sensitivity",
      "kmer_length",
      "max_self_map",
      "alignment_kernel",
      "extra_band_width",
      "max_alignment_attempts",
    ]);
    assert_eq!(args_under(&PangraphBuildArgs::command(), "Alignment"), alignment);
    assert_eq!(args_under(&PangraphMergeArgs::command(), "Alignment"), alignment);
  }
}
