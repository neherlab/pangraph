use crate::commands::root_args::PangraphArgs;
use eyre::{Report, WrapErr};
use regex::Regex;
use std::borrow::Cow;

pub fn print_help_markdown() -> Result<(), Report> {
  let help = clap_markdown::help_markdown::<PangraphArgs>();

  let help = replace(&help, "# Command-Line Help for `pangraph`", "")?;

  let help = replace(
    &help,
    "This document contains the help content for the `pangraph` command-line program.",
    r#"
This document contains the automatically generated reference documentation for command-line arguments of the latest version of Pangraph CLI.

If you have Pangraph CLI installed, you can type `pangraph --help` to read the latest documentation for your installed version of Pangraph. To generate this document in markdown format, run `pangraph help-markdown > reference.md`
  "#,
  )?;

  let help = replace(&help, "(.*)— REMOVED(.*)", "")?;
  let help = replace(&help, "(.*)— RENAMED(.*)", "")?;

  let help = replace(
    &help,
    r#"(?<orig>\* `--server <SERVER>` — Use custom dataset server)\n\n  Default value: .*"#,
    "$orig",
  )?;

  let help = replace(
    &help,
    r#"(?<orig>\* `-j`, `--jobs <JOBS>` — Number of processing jobs. If not specified, all available CPU threads will be used)\n\n  Default value: .*"#,
    "$orig",
  )?;

  let help = replace(&help, "", "")?;

  println!("{help}");
  Ok(())
}

fn replace<'t>(text: &'t str, what: &str, with_what: &str) -> Result<Cow<'t, str>, Report> {
  let res = Regex::new(what)
    .wrap_err_with(|| format!("When compiling regex: {what}"))?
    .replace_all(text, with_what);
  Ok(res)
}
