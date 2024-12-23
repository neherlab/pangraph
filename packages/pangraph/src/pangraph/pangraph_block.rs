use crate::io::seq::reverse_complement;
use crate::pangraph::edits::Edit;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_node::NodeId;
use crate::pangraph::pangraph_path::PathId;
use crate::utils::collections::has_duplicates;
use derive_more::{Display, From};
use eyre::{Report, WrapErr};
use getset::{CopyGetters, Getters};
use maplit::btreemap;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::hash::Hash;

#[derive(
  Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize, JsonSchema,
)]
pub struct BlockId(pub usize);

impl BlockId {
  pub fn from_str(block_id_str: &String) -> Result<Self, Report> {
    let block_id = block_id_str
      .parse::<usize>()
      .wrap_err_with(|| format!("When parsing Block ID: expected unsigned integer, but got '{block_id_str}'"))?;
    Ok(Self(block_id))
  }
}

#[must_use]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, Getters, CopyGetters, JsonSchema)]
pub struct PangraphBlock {
  #[get_copy = "pub"]
  id: BlockId,
  consensus: String,
  alignments: BTreeMap<NodeId, Edit>,
}

impl PangraphBlock {
  pub fn from_consensus(consensus: impl Into<String>, block_id: BlockId, nid: NodeId) -> Self {
    PangraphBlock::new(block_id, consensus, btreemap! {nid => Edit::empty()})
  }

  pub fn new(block_id: BlockId, consensus: impl Into<String>, alignments: BTreeMap<NodeId, Edit>) -> Self {
    let consensus = consensus.into();
    let id = block_id;
    Self {
      id,
      consensus,
      alignments,
    }
  }

  pub fn reverse_complement(&self) -> Result<Self, Report> {
    let rev_cons = reverse_complement(&self.consensus)?;

    let len = self.consensus_len();

    let rev_aln = self
      .alignments()
      .iter()
      .map(|(&nid, e)| Ok((nid, e.reverse_complement(len)?)))
      .collect::<Result<_, Report>>()?;

    Ok(Self::new(self.id(), rev_cons, rev_aln))
  }

  pub fn depth(&self) -> usize {
    self.alignments.len()
  }

  pub fn consensus(&self) -> &str {
    self.consensus.as_str()
  }

  pub fn consensus_mut(&mut self) -> &mut str {
    self.consensus.as_mut_str()
  }

  pub fn set_consensus(&mut self, consensus: impl Into<String>) {
    self.consensus = consensus.into();
  }

  #[allow(unsafe_code)]
  pub fn set_consensus_char(&mut self, pos: usize, c: char) {
    debug_assert!(self.consensus.is_ascii());
    // FIXME: we REALLY should not use strings for sequences anymore
    // SAFETY: safe as long as we have ASCII strings
    let target = unsafe { &mut self.consensus.as_bytes_mut()[pos] };
    *target = c as u8;
  }

  pub fn consensus_len(&self) -> usize {
    self.consensus.len()
  }

  pub fn alignment(&self, nid: NodeId) -> &Edit {
    &self.alignments[&nid]
  }

  pub fn alignment_mut(&mut self, nid: NodeId) -> &mut Edit {
    self.alignments.get_mut(&nid).unwrap()
  }

  pub fn alignment_keys(&self) -> BTreeSet<NodeId> {
    self.alignments.keys().copied().collect()
  }

  pub fn alignments(&self) -> &BTreeMap<NodeId, Edit> {
    &self.alignments
  }

  pub fn alignments_mut(&mut self) -> &mut BTreeMap<NodeId, Edit> {
    &mut self.alignments
  }

  pub fn alignment_insert(&mut self, node_id: NodeId, edit: Edit) -> Option<Edit> {
    self.alignments.insert(node_id, edit)
  }

  pub fn alignment_remove(&mut self, node_id: NodeId) -> Edit {
    self.alignments.remove(&node_id).unwrap()
  }

  pub fn isolates<'a>(&'a self, graph: &'a Pangraph) -> impl Iterator<Item = PathId> + 'a {
    self.alignments.keys().map(|&node_id| graph.nodes[&node_id].path_id())
  }

  pub fn is_duplicated(&self, graph: &Pangraph) -> bool {
    has_duplicates(self.isolates(graph))
  }
}
