use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_node::PangraphNode;
use crate::pangraph::strand::Strand;
use crate::utils::id::Id;
use derive_more::{Display, From};
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct PathId(pub usize);

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct PangraphPath {
  pub name: String,
  pub nodes: Vec<PangraphNode>,
  pub tot_len: usize,
  pub circular: bool,
}

impl Id<PathId> for PangraphPath {}

impl PangraphPath {
  pub fn new(name: impl AsRef<str>, block_id: BlockId, strand: Strand, circular: bool) -> Self {
    let path_id = PathId(0); // FIXME: self.id();
    Self {
      name: name.as_ref().to_owned(),
      nodes: vec![PangraphNode::new(block_id, path_id, strand, (0, 0))],
      tot_len: 0,
      circular,
    }
  }
}
