use crate::io::fasta::FastaRecord;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_parse;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::pangraph::strand::Strand;
use crate::utils::map_merge::{map_merge, ConflictResolution};
use eyre::{Report, WrapErr};
use maplit::btreemap;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::str::FromStr;

#[derive(Clone, Debug, Default, Serialize, Deserialize, Hash)]
pub struct Pangraph {
  pub paths: BTreeMap<PathId, PangraphPath>,
  pub blocks: BTreeMap<BlockId, PangraphBlock>,
  pub nodes: BTreeMap<NodeId, PangraphNode>,
}

impl Pangraph {
  pub fn singleton(fasta: FastaRecord, strand: Strand, circular: bool) -> Self {
    let tot_len = fasta.seq.len();
    let node_id = NodeId(fasta.index);
    let block = PangraphBlock::from_consensus(fasta.seq, node_id);
    let path_id = PathId(fasta.index);
    let node = PangraphNode::new(Some(node_id), block.id(), path_id, strand, (0, 0));
    let path = PangraphPath::new(Some(path_id), [node.id()], tot_len, circular, Some(fasta.seq_name));
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
    self.blocks.values().map(|block| block.consensus())
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

    self.blocks = map_merge(
      &self.blocks,
      &u.b_new.iter().map(|b| (b.id(), b.clone())).collect(),
      ConflictResolution::Custom(|(kl, vl), (kr, vr)| panic!("Conflicting key: '{kl}'")),
    );

    for (old_node_id, new_nodes) in &u.n_new {
      let path_id = self.nodes[&old_node_id].path_id();

      let path = self.paths.get_mut(&path_id).unwrap(); // FIXME

      let old_idx = path.nodes.iter().position(|node_id| node_id == old_node_id).unwrap();

      path.nodes.remove(old_idx);

      let new_ids: Vec<NodeId> = new_nodes.iter().map(|n| n.id()).collect();
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
  use crate::pangraph::edits::Edit;
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::pangraph_path::PangraphPath;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
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
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(1), PathId(1), Forward,  (0, 0)), // FIXME
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), BlockId(1), PathId(3), Forward,  (0, 0)), // FIXME
      NodeId(3) => PangraphNode::new(Some(NodeId(3)), BlockId(2), PathId(1), Forward,  (0, 0)), // FIXME
      NodeId(4) => PangraphNode::new(Some(NodeId(4)), BlockId(2), PathId(2), Forward,  (0, 0)), // FIXME
      NodeId(5) => PangraphNode::new(Some(NodeId(5)), BlockId(2), PathId(3), Reverse, (0, 0)), // FIXME
      NodeId(6) => PangraphNode::new(Some(NodeId(6)), BlockId(3), PathId(1), Forward,  (0, 0)), // FIXME
      NodeId(7) => PangraphNode::new(Some(NodeId(7)), BlockId(3), PathId(2), Forward,  (0, 0)), // FIXME
      NodeId(8) => PangraphNode::new(Some(NodeId(8)), BlockId(3), PathId(3), Forward,  (0, 0)) // FIXME
    };

    let blocks = btreemap! {
      BlockId(1) => PangraphBlock::new(Some(BlockId(1)), "1",
        btreemap!{ NodeId(1) => Edit::empty(), NodeId(2) => Edit::empty() }),
      BlockId(2) => PangraphBlock::new(Some(BlockId(2)), "2",
        btreemap!{ NodeId(3) => Edit::empty(), NodeId(4) => Edit::empty(), NodeId(5) => Edit::empty() }),
      BlockId(3) => PangraphBlock::new(Some(BlockId(3)), "3",
        btreemap!{ NodeId(6) => Edit::empty(), NodeId(7) => Edit::empty(), NodeId(8) => Edit::empty(), }),
    };

    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1), NodeId(3), NodeId(6)], 0, false, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [NodeId(4), NodeId(7)           ], 0, false, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [NodeId(2), NodeId(5), NodeId(8)], 0, false, None),
    };

    let mut G = Pangraph {
      paths: paths.clone(),
      blocks: blocks.clone(),
      nodes: nodes.clone(),
    };

    let new_nodes = btreemap! {
      NodeId(9)  => PangraphNode::new(Some(NodeId(9)),  BlockId(4), PathId(1), Forward,  (0, 0)),
      NodeId(10) => PangraphNode::new(Some(NodeId(10)), BlockId(5), PathId(1), Reverse, (0, 0)),
      NodeId(11) => PangraphNode::new(Some(NodeId(11)), BlockId(4), PathId(2), Forward,  (0, 0)),
      NodeId(12) => PangraphNode::new(Some(NodeId(12)), BlockId(5), PathId(2), Reverse, (0, 0)),
      NodeId(13) => PangraphNode::new(Some(NodeId(13)), BlockId(4), PathId(3), Reverse, (0, 0)),
      NodeId(14) => PangraphNode::new(Some(NodeId(14)), BlockId(5), PathId(3), Forward,  (0, 0)),
    };

    let new_blocks = btreemap! {
      BlockId(4) => PangraphBlock::new(Some(BlockId(4)), "4", btreemap!{}),
      BlockId(5) => PangraphBlock::new(Some(BlockId(5)), "5", btreemap!{}),
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
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1),  NodeId(9),  NodeId(10),  NodeId(6)], 0, false, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [NodeId(11), NodeId(12), NodeId(7)             ], 0, false, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [NodeId(2),  NodeId(14), NodeId(13),  NodeId(8)], 0, false, None),
    };
    assert_eq!(G.paths, expected_paths);

    let expected_nodes = btreemap! {
      NodeId(1) => nodes[&NodeId(1)].clone(),
      NodeId(2) => nodes[&NodeId(2)].clone(),
      NodeId(6) => nodes[&NodeId(6)].clone(),
      NodeId(7) => nodes[&NodeId(7)].clone(),
      NodeId(8) => nodes[&NodeId(8)].clone(),
      NodeId(9) => new_nodes[&NodeId(9)].clone(),
      NodeId(10) => new_nodes[&NodeId(10)].clone(),
      NodeId(11) => new_nodes[&NodeId(11)].clone(),
      NodeId(12) => new_nodes[&NodeId(12)].clone(),
      NodeId(13) => new_nodes[&NodeId(13)].clone(),
      NodeId(14) => new_nodes[&NodeId(14)].clone(),
    };
    assert_eq!(G.nodes, expected_nodes);
  }
}
