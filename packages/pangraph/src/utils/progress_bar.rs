use crate::utils::global_init::INDICATIF;
use atty::{is as is_tty, Stream};
use eyre::Report;
use indicatif::{ProgressBar as ProgressBarBase, ProgressStyle};
use std::borrow::Cow;
use std::time::Duration;

pub struct ProgressBar {
  pb: Option<ProgressBarBase>,
}

impl ProgressBar {
  pub fn new(n_total: usize) -> Result<Self, Report> {
    let pb = is_tty(Stream::Stdout).then(|| {
      let pb = ProgressBarBase::new(n_total as u64);
      pb.enable_steady_tick(Duration::from_secs(1));
      pb.set_style(
        ProgressStyle::with_template(
          "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{human_len} mergers ({eta})",
        )
        .unwrap()
        .progress_chars("#>-"),
      );
      let pb = INDICATIF.get().unwrap().add(pb);
      pb
    });
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
