use crate::make_error;
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use std::ffi::OsStr;
use std::fmt::Debug;
use std::process::{Command, Stdio};

pub fn subprocess_with_args<I, S>(command: S, args: I) -> Result<String, Report>
where
  I: IntoIterator<Item = S>,
  S: AsRef<OsStr> + Debug,
{
  subprocess(&command)
    .wrap_err_with(|| format!("When running command {command:?}"))
    .with_section(|| {
      format!("Make sure you have {command:#?} installed and available in system $PATH.").header("Suggestion: ")
    })?;

  let child = Command::new(&command)
    .args(args)
    .stdout(Stdio::piped())
    .stderr(Stdio::piped())
    .spawn()?;

  let output = child.wait_with_output()?;

  if !output.status.success() {
    return make_error!(
      "Subprocess failed:\n  Exit code: '{}'\n  stderr: {}\n  stdout: {}",
      output.status.code().unwrap_or(0),
      String::from_utf8(output.stderr)?,
      String::from_utf8(output.stdout)?
    );
  }

  Ok(String::from_utf8(output.stdout)?)
}

pub fn subprocess<S>(command: S) -> Result<String, Report>
where
  S: AsRef<OsStr>,
{
  let child = Command::new(command)
    .stdout(Stdio::piped())
    .stderr(Stdio::piped())
    .spawn()?;

  let output = child.wait_with_output()?;

  if !output.status.success() {
    return make_error!(
      "Subprocess failed:\n  Exit code: '{}'\n  stderr: {}\n  stdout: {}",
      output.status.code().unwrap_or(0),
      String::from_utf8(output.stderr)?,
      String::from_utf8(output.stdout)?
    );
  }

  Ok(String::from_utf8(output.stdout)?)
}

/// Creates an argument for Command::arg from name and value.
/// Conditionally added argument need to be a vec of 0, 1 or 2 strings.
pub fn create_arg_optional(name: impl AsRef<str>, value: &Option<impl ToString>) -> Vec<String> {
  value
    .as_ref()
    .map(|kmer_len| vec![name.as_ref().to_owned(), kmer_len.to_string()])
    .unwrap_or_default()
}
