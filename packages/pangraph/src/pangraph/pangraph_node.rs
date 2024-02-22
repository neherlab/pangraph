use crate::pangraph::strand::Strand;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PangraphNode {
  pub id: usize,
  pub block_id: usize,
  pub path_id: usize,
  pub strand: Strand,
  pub position: (usize, usize),
}
