use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, SmartDefault)]
pub enum Strand {
  #[default]
  #[serde(rename = "+")]
  Forward,
  #[serde(rename = "-")]
  Reverse,
}
