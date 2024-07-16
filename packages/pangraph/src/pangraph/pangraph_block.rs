use crate::align::map_variations::map_variations;
use crate::pangraph::edits::Edit;
use crate::pangraph::pangraph_node::NodeId;
use crate::utils::id::id;
use derive_more::{Display, From};
use eyre::{Report, WrapErr};
use getset::{CopyGetters, Getters};
use maplit::btreemap;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::hash::Hash;

#[derive(Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct BlockId(pub usize);

impl BlockId {
  pub fn from_str(block_id_str: &String) -> Result<Self, Report> {
    let block_id = block_id_str
      .parse::<usize>()
      .wrap_err_with(|| format!("When parsing Block ID: expected unsigned integer, but got '{block_id_str}'"))?;
    Ok(Self(block_id))
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, Getters, CopyGetters)]
pub struct PangraphBlock {
  #[get_copy = "pub"]
  id: BlockId,
  consensus: String,
  alignments: BTreeMap<NodeId, Edit>,
}

impl PangraphBlock {
  pub fn from_consensus(consensus: impl Into<String>, nid: NodeId) -> Self {
    PangraphBlock::new(None, consensus, btreemap! {nid => Edit::empty()})
  }

  pub fn new(block_id: Option<BlockId>, consensus: impl Into<String>, alignments: BTreeMap<NodeId, Edit>) -> Self {
    let consensus = consensus.into();
    let id = block_id.unwrap_or_else(|| id((&consensus, &alignments)));
    Self {
      id,
      consensus,
      alignments,
    }
  }

  pub fn refresh_id(&mut self) {
    self.id = id((&self.consensus, &self.alignments));
  }

  pub fn depth(&self) -> usize {
    self.alignments.len()
  }

  pub fn consensus(&self) -> &str {
    self.consensus.as_str()
  }

  pub fn consensus_len(&self) -> usize {
    self.consensus.len()
  }

  pub fn alignment(&self, nid: NodeId) -> &Edit {
    &self.alignments[&nid]
  }

  pub fn alignment_keys(&self) -> BTreeSet<NodeId> {
    self.alignments.keys().copied().collect()
  }

  pub fn alignments(&self) -> &BTreeMap<NodeId, Edit> {
    &self.alignments
  }

  pub fn alignment_insert(&mut self, node_id: NodeId, edit: Edit) -> Option<Edit> {
    self.alignments.insert(node_id, edit)
  }
}
