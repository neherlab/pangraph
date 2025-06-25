use crate::utils::global_init::INDICATIF;
use atty::{Stream, is as is_tty};
use eyre::Report;
use indicatif::{ProgressBar as ProgressBarBase, ProgressStyle};
use std::borrow::Cow;
use std::time::Duration;

pub struct ProgressBar {
  pb: Option<ProgressBarBase>,
}

impl ProgressBar {
  pub fn new(n_total: usize, deactivate: bool) -> Result<Self, Report> {
    if deactivate || (n_total <= 1) {
      return Ok(Self { pb: None });
    }
    let pb = if is_tty(Stream::Stdout) {
      let pb = ProgressBarBase::new(n_total as u64);
      pb.enable_steady_tick(Duration::from_secs(1));
      let style = ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {human_pos}/{human_len} mergers",
      )
      .map_err(|e| eyre::eyre!("Failed to create progress bar template: {}", e))?
      .progress_chars("#>-");
      pb.set_style(style);
      let indicatif = INDICATIF.get()
        .ok_or_else(|| eyre::eyre!("Failed to get INDICATIF instance"))?;
      Some(indicatif.add(pb))
    } else {
      None
    };
    Ok(Self { pb })
  }

  pub fn inc(&self, delta: u64) {
    if let Some(pb) = &self.pb {
      pb.inc(delta);
    }
  }

  pub fn finish_with_message(&self, msg: impl Into<Cow<'static, str>>) {
    if let Some(pb) = &self.pb {
      pb.finish_with_message(msg);
    }
  }
}
