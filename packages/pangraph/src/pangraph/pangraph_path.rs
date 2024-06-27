use crate::pangraph::pangraph_node::NodeId;
use crate::utils::id::Id;
use derive_more::{Display, From};
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct PathId(pub usize);

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct PangraphPath {
  // pub name: BlockId,
  pub nodes: Vec<NodeId>,
  pub tot_len: usize,
  pub circular: bool,
}

impl Id<PathId> for PangraphPath {}

impl PangraphPath {
  pub fn new(/* name: BlockId, */ nodes: impl Into<Vec<NodeId>>, tot_len: usize, circular: bool) -> Self {
    Self {
      // name,
      nodes: nodes.into(),
      tot_len,
      circular,
    }
  }
}
