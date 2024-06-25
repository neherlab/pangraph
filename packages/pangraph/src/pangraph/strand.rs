use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, SmartDefault)]
pub enum Strand {
  #[default]
  #[serde(rename = "+")]
  Forward,
  #[serde(rename = "-")]
  Reverse,
}
