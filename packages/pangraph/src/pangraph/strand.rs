use crate::{make_error, make_report};
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[must_use]
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, SmartDefault)]
pub enum Strand {
  #[default]
  #[serde(rename = "+")]
  Forward,
  #[serde(rename = "-")]
  Reverse,
}

impl Strand {
  pub fn is_forward(&self) -> bool {
    matches!(self, Strand::Forward)
  }

  pub fn is_reverse(&self) -> bool {
    matches!(self, Strand::Reverse)
  }

  pub fn reverse(&self) -> Strand {
    match self {
      Strand::Forward => Strand::Reverse,
      Strand::Reverse => Strand::Forward,
    }
  }

  pub fn from_char(c: char) -> Result<Self, Report> {
    match c {
      '+' => Ok(Strand::Forward),
      '-' => Ok(Strand::Reverse),
      _ => make_error!("Invalid strand character: {c}"),
    }
  }

  pub fn from_str(s: &str) -> Result<Self, Report> {
    if s.len() != 1 {
      make_error!("Invalid strand string: '{s}'")
    } else {
      let c = s.chars().next().ok_or_else(|| make_report!("Empty strand string"))?;
      Strand::from_char(c)
    }
  }
}
