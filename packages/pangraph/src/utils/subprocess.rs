use crate::make_error;
use eyre::Report;
use std::ffi::OsStr;
use std::process::{Command, Stdio};

pub fn subprocess<I, S>(command: S, args: I) -> Result<String, Report>
where
  I: IntoIterator<Item = S>,
  S: AsRef<OsStr>,
{
  let child = Command::new(command)
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
