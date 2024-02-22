use crate::align::alignment::Alignment;
use crate::utils::id::random_id;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct PangraphBlock {
  pub id: usize,
  pub consensus: String,
  pub alignment: Alignment,
}

impl PangraphBlock {
  pub fn new(consensus: String) -> Self {
    Self {
      id: random_id(),
      consensus,
      alignment: Alignment::default(),
    }
  }
}
