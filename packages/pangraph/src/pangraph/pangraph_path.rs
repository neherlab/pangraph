use crate::pangraph::pangraph_node::NodeId;
use crate::utils::id::Id;
use derive_more::{Display, From};
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct PathId(pub usize);

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct PangraphPath {
  pub name: String,
  pub nodes: Vec<NodeId>,
  pub tot_len: usize,
  pub circular: bool,
}

impl Id<PathId> for PangraphPath {}

impl PangraphPath {
  pub fn new(name: impl AsRef<str>, nodes: &[NodeId], circular: bool) -> Self {
    Self {
      name: name.as_ref().to_owned(),
      nodes: nodes.to_owned(),
      tot_len: 0,
      circular,
    }
  }
}
