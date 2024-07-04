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
}

#[must_use]
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, SmartDefault)]
pub enum Orientation {
  #[default]
  #[serde(rename = "+")]
  Forward,
  #[serde(rename = "-")]
  Reverse,
}

impl Orientation {
  pub fn is_forward(&self) -> bool {
    matches!(self, Orientation::Forward)
  }

  pub fn is_reverse(&self) -> bool {
    matches!(self, Orientation::Reverse)
  }

  pub fn reverse(&self) -> Orientation {
    match self {
      Orientation::Forward => Orientation::Reverse,
      Orientation::Reverse => Orientation::Forward,
    }
  }
}
