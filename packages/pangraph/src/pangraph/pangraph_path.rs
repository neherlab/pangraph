use crate::pangraph::pangraph_node::PangraphNode;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PangraphPath {
  pub id: usize,
  pub nodes: Vec<PangraphNode>,
  pub tot_len: usize, // TODO: potentially can be calculated from other fields?
  // pub offset: isize,
  pub circular: bool,
  // pub positions: BTreeMap<String, usize>,
}
