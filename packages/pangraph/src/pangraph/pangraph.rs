use crate::io::fasta::FastaRecord;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_parse;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::utils::id::Id;
use eyre::{Report, WrapErr};
use maplit::btreemap;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::str::FromStr;

#[derive(Clone, Serialize, Deserialize, Debug, Hash)]
pub struct Pangraph {
  pub paths: BTreeMap<PathId, PangraphPath>,
  pub blocks: BTreeMap<BlockId, PangraphBlock>,
  pub nodes: BTreeMap<NodeId, PangraphNode>,
}

impl Pangraph {
  pub fn singleton(fasta: FastaRecord, strand: bool, circular: bool) -> Self {
    let block = PangraphBlock::from_consensus(fasta.seq);

    // FIXME: Paths and Nodes depend on ids of each other - chicken and egg problem. How to have simultaneously:
    //
    // - correct hash ids (without ad-hoc external id calculation and partially initialized objects)
    //   AND
    // - paths and nodes cross-reference each other
    //
    // Use pointers instead?
    let path_id = PathId(0);
    let node = PangraphNode::new(block.id(), path_id, strand, (0, 0));

    let tot_len = 0; // FIXME
    let path = PangraphPath::new(/*fasta.seq_name, */ [node.id()], tot_len, circular);

    Self {
      blocks: btreemap! {block.id() => block},
      paths: btreemap! {path.id() => path},
      nodes: btreemap! {node.id() => node},
    }
  }

  pub fn from_path<P: AsRef<Path>>(filepath: &Option<P>) -> Result<Self, Report> {
    let reader = open_file_or_stdin(filepath)?;
    let data = read_reader_to_string(reader).wrap_err("When reading Pangraph JSON")?;
    Self::from_str(&data).wrap_err("When parsing Pangraph JSON")
  }

  pub fn to_string_pretty(&self) -> Result<String, Report> {
    let mut tree_str = serde_json::to_string_pretty(self)?;
    tree_str += "\n";
    Ok(tree_str)
  }

  pub fn consensuses(&self) -> impl Iterator<Item = &str> {
    self.blocks.values().map(|block| block.consensus.as_str())
  }

  pub fn update(&mut self, u: &GraphUpdate) {
    // Consistency check: node ids
    let old_nodes_set_from_graph: BTreeSet<NodeId> = self.blocks[&u.b_old_id].alignment_keys();
    let old_nodes_set_from_update: BTreeSet<NodeId> = u.n_new.keys().copied().collect();

    assert_eq!(
      old_nodes_set_from_graph, old_nodes_set_from_update,
      "old nodes mismatch: {old_nodes_set_from_graph:#?} != {old_nodes_set_from_update:#?}",
    );

    self.blocks.remove(&u.b_old_id);

    for b in &u.b_new {
      self.blocks.insert(b.id(), b.clone());
    }

    for (old_node_id, new_nodes) in &u.n_new {
      let path_id = self.nodes[&old_node_id].path_id();

      let path = self.paths.get_mut(&path_id).unwrap(); // FIXME

      let old_idx = path.nodes.iter().position(|node_id| node_id == old_node_id).unwrap();

      path.nodes.remove(old_idx);

      let new_ids: Vec<NodeId> = new_nodes.iter().map(Id::id).collect();
      path.nodes.splice(old_idx..old_idx, new_ids);

      self.nodes.remove(old_node_id);

      for n in new_nodes {
        self.nodes.insert(n.id(), n.to_owned());
      }
    }
  }
}

impl FromStr for Pangraph {
  type Err = Report;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    json_parse(s).wrap_err("When parsing Pangraph JSON contents")
  }
}

#[derive(Debug)]
pub struct GraphUpdate {
  pub b_old_id: BlockId,
  pub b_new: Vec<PangraphBlock>,
  pub n_new: BTreeMap<NodeId, Vec<PangraphNode>>,
  // nb: node list is already in the order of the new path
}

#[cfg(test)]
mod tests {
  #![allow(non_snake_case, clippy::redundant_clone)]

  use super::*;
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::pangraph_path::PangraphPath;
  use maplit::btreemap;

  #[test]
  fn test_graph_update() {
    // graph
    // p1 -> [b1+,b2+,b3+]
    // p2 -> [b2+,b3+]
    // p3 -> [b1+,b2-,b3+]
    // update
    // b2+ -> [b4+, b5-]

    let nodes = btreemap! {
      NodeId(1) => PangraphNode::new(/* 1, */ BlockId(1), PathId(1), true, (0, 0)), // FIXME
      NodeId(2) => PangraphNode::new(/* 2, */ BlockId(1), PathId(3), true, (0, 0)), // FIXME
      NodeId(3) => PangraphNode::new(/* 3, */ BlockId(2), PathId(1), true, (0, 0)), // FIXME
      NodeId(4) => PangraphNode::new(/* 4, */ BlockId(2), PathId(2), true, (0, 0)), // FIXME
      NodeId(5) => PangraphNode::new(/* 5, */ BlockId(2), PathId(3), false, (0, 0)), // FIXME
      NodeId(6) => PangraphNode::new(/* 6, */ BlockId(3), PathId(1), true, (0, 0)), // FIXME
      NodeId(7) => PangraphNode::new(/* 7, */ BlockId(3), PathId(2), true, (0, 0)), // FIXME
      NodeId(8) => PangraphNode::new(/* 8, */ BlockId(3), PathId(3), true, (0, 0)) // FIXME
    };

    let blocks = btreemap! {
      BlockId(1) => PangraphBlock::new(/* 1, */ "1", btreemap!{}), // {1: None, 2: None}
      BlockId(2) => PangraphBlock::new(/* 2, */ "2", btreemap!{}), // {3: None, 4: None, 5: None}
      BlockId(3) => PangraphBlock::new(/* 3, */ "3", btreemap!{}), // {6: None, 7: None, 8: None}
    };

    let paths = btreemap! {
      PathId(1) => PangraphPath::new(/* 1, */ [NodeId(1), NodeId(3), NodeId(6)], 0, false),
      PathId(2) => PangraphPath::new(/* 2, */ [NodeId(4), NodeId(7)], 0, false),
      PathId(3) => PangraphPath::new(/* 3, */ [NodeId(2), NodeId(5), NodeId(8)], 0, false)
    };

    let mut G = Pangraph {
      paths: paths.clone(),
      blocks: blocks.clone(),
      nodes: nodes.clone(),
    };

    let new_nodes = btreemap! {
      NodeId(9)  => PangraphNode::new(/* 9,  */ BlockId(4), PathId(1), true,  (0, 0)),
      NodeId(10) => PangraphNode::new(/* 10, */ BlockId(5), PathId(1), false, (0, 0)),
      NodeId(11) => PangraphNode::new(/* 11, */ BlockId(4), PathId(2), true,  (0, 0)),
      NodeId(12) => PangraphNode::new(/* 12, */ BlockId(5), PathId(2), false, (0, 0)),
      NodeId(13) => PangraphNode::new(/* 13, */ BlockId(4), PathId(3), false, (0, 0)),
      NodeId(14) => PangraphNode::new(/* 14, */ BlockId(5), PathId(3), true,  (0, 0))
    };

    let new_blocks = btreemap! {
      BlockId(4) => PangraphBlock::new(/* 4, */ "4", btreemap!{}),
      BlockId(5) => PangraphBlock::new(/* 5, */ "5", btreemap!{})
    };

    let update = GraphUpdate {
      b_old_id: BlockId(2),
      b_new: vec![new_blocks[&BlockId(4)].clone(), new_blocks[&BlockId(5)].clone()],
      n_new: btreemap! {
        NodeId(3) => vec![new_nodes[&NodeId(9) ].clone(), new_nodes[&NodeId(10)].clone()],
        NodeId(4) => vec![new_nodes[&NodeId(11)].clone(), new_nodes[&NodeId(12)].clone()],
        NodeId(5) => vec![new_nodes[&NodeId(14)].clone(), new_nodes[&NodeId(13)].clone()],
      },
    };

    G.update(&update);

    let expected_blocks = btreemap! {
      BlockId(1) => blocks[&BlockId(1)].clone(),
      BlockId(3) => blocks[&BlockId(3)].clone(),
      BlockId(4) => new_blocks[&BlockId(4)].clone(),
      BlockId(5) => new_blocks[&BlockId(5)].clone(),
    };
    assert_eq!(G.blocks, expected_blocks);

    let expected_paths = btreemap! {
      PathId(1) => PangraphPath::new(/* 1, */ [NodeId(1),  NodeId(9),  NodeId(10),  NodeId(6)], 0, false),
      PathId(2) => PangraphPath::new(/* 2, */ [NodeId(11), NodeId(12), NodeId(7)             ], 0, false),
      PathId(3) => PangraphPath::new(/* 3, */ [NodeId(2),  NodeId(14), NodeId(13),  NodeId(8)], 0, false),
    };
    assert_eq!(G.paths, expected_paths);

    let expected_nodes = btreemap! {
      NodeId(1) => nodes[&NodeId(1)].clone(),
      NodeId(2) => nodes[&NodeId(2)].clone(),
      NodeId(6) => nodes[&NodeId(6)].clone(),
      NodeId(7) => nodes[&NodeId(7)].clone(),
      NodeId(8) => nodes[&NodeId(8)].clone(),
      NodeId(9) => nodes[&NodeId(9)].clone(),
      NodeId(10) => nodes[&NodeId(10)].clone(),
      NodeId(11) => nodes[&NodeId(11)].clone(),
      NodeId(12) => nodes[&NodeId(12)].clone(),
      NodeId(13) => nodes[&NodeId(13)].clone(),
      NodeId(14) => nodes[&NodeId(14)].clone(),
    };
    assert_eq!(G.nodes, expected_nodes);
  }
}
