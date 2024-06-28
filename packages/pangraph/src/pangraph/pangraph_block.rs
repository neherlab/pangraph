use crate::align::map_variations::map_variations;
use crate::pangraph::edits::Edit;
use crate::pangraph::pangraph_node::NodeId;
use crate::utils::id::id;
use derive_more::{Display, From};
use eyre::Report;
use getset::{CopyGetters, Getters};
use maplit::btreemap;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::hash::Hash;

#[derive(Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct BlockId(pub usize);

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, Getters, CopyGetters)]
pub struct PangraphBlock {
  #[get_copy = "pub"]
  id: BlockId,
  consensus: String,
  alignments: BTreeMap<NodeId, Edit>,
}

impl PangraphBlock {
  pub fn from_consensus(consensus: impl Into<String>) -> Self {
    PangraphBlock::new(None, consensus, btreemap! {})
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

  /// Append a sequence to a block, and assign the given node_id
  pub fn append_sequence(&mut self, new_seq: impl AsRef<str>, new_node_id: NodeId) -> Result<(), Report> {
    // Calculate the set of edits
    let edits = map_variations(&self.consensus, new_seq)?;
    // Append the edits to the block alignment dictionary
    self.alignments.insert(new_node_id, edits);
    self.refresh_id();
    Ok(())
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
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::pangraph::edits::{Del, Ins, Sub};
  use pretty_assertions::assert_eq;

  #[test]
  fn test_block_append_sequence_simple_case() {
    let consensus = o!("ACTTTGCGTATTTACTATA");

    let E_0 = Edit {
      dels: vec![Del::new(2, 3)],
      subs: vec![],
      inss: vec![],
    };
    let E_1 = Edit {
      inss: vec![Ins::new(16, "CCC")],
      subs: vec![],
      dels: vec![],
    };
    let E_2 = Edit {
      subs: vec![Sub::new(11, 'C')],
      dels: vec![],
      inss: vec![],
    };

    let alignments = BTreeMap::from([(NodeId(0), E_0), (NodeId(1), E_1), (NodeId(2), E_2)]);

    let mut block = PangraphBlock::new(None, consensus, alignments);

    let new_seq = o!("ACTTTGCGGATTTACTATA");
    block.append_sequence(new_seq, NodeId(0)).unwrap();

    let E_new_expected = Edit {
      subs: vec![Sub::new(8, 'G')],
      dels: vec![],
      inss: vec![],
    };

    assert_eq!(block.alignments[&NodeId(0)], E_new_expected);
  }
}
