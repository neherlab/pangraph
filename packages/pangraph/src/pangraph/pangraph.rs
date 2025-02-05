use crate::io::fasta::FastaRecord;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_read_str;
use crate::make_internal_report;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::pangraph::strand::Strand;
use crate::representation::seq::Seq;
use crate::utils::map_merge::{map_merge, ConflictResolution};
use eyre::{Report, WrapErr};
use maplit::btreemap;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::str::FromStr;

#[derive(Clone, Debug, Default, Serialize, Deserialize, Hash, PartialEq, Eq, JsonSchema)]
pub struct Pangraph {
  pub paths: BTreeMap<PathId, PangraphPath>,
  pub blocks: BTreeMap<BlockId, PangraphBlock>,
  pub nodes: BTreeMap<NodeId, PangraphNode>,
}

impl Pangraph {
  pub fn singleton(fasta: FastaRecord, strand: Strand, circular: bool) -> Self {
    let tot_len = fasta.seq.len();
    let node_id = NodeId(fasta.index);
    let block_id = BlockId(fasta.index);
    let block = PangraphBlock::from_consensus(fasta.seq, block_id, node_id);
    let path_id = PathId(fasta.index);
    let node_position = if circular { (0, 0) } else { (0, tot_len) }; // path wraps around if circular
    let node = PangraphNode::new(Some(node_id), block.id(), path_id, strand, node_position);
    let path = PangraphPath::new(
      Some(path_id),
      [node.id()],
      tot_len,
      circular,
      Some(fasta.seq_name),
      fasta.desc,
    );
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

  pub fn consensuses(&self) -> impl Iterator<Item = &Seq> {
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
      let path_id = self.nodes[old_node_id].path_id();

      let path = self.paths.get_mut(&path_id).unwrap(); // FIXME

      let old_idx = path.nodes.iter().position(|node_id| node_id == old_node_id).unwrap();

      // only one such nodes
      debug_assert!(path.nodes.iter().filter(|node_id| *node_id == old_node_id).count() == 1);

      path.nodes.remove(old_idx);

      let new_ids: Vec<NodeId> = new_nodes.iter().map(|n| n.id()).collect();
      path.nodes.splice(old_idx..old_idx, new_ids);

      self.nodes.remove(old_node_id);

      for n in new_nodes {
        self.nodes.insert(n.id(), n.to_owned());
      }
    }
  }

  #[allow(unused_must_use)]
  pub fn remove_path(&mut self, pid: PathId) {
    if let Some(path) = self.paths.remove(&pid) {
      for nid in path.nodes {
        if let Some(node) = self.nodes.remove(&nid) {
          if let Some(block) = self.blocks.get_mut(&node.block_id()) {
            block.alignment_remove(nid);
          }
        }
      }
    }

    // remove empty blocks
    let empty_blocks: Vec<BlockId> = self
      .blocks
      .iter()
      .filter(|(bid, block)| block.alignments().is_empty())
      .map(|(bid, _)| *bid)
      .collect();

    for bid in empty_blocks {
      self.blocks.remove(&bid);
    }
  }

  #[cfg(any(test, debug_assertions))]
  pub fn sanity_check(&self) -> Result<(), Report> {
    for (node_id, node) in &self.nodes {
      if !self.blocks.contains_key(&node.block_id()) {
        return Err(eyre::eyre!("Block {} not found in graph", node.block_id()));
      }
      let block = &self.blocks[&node.block_id()];

      if !self.paths.contains_key(&node.path_id()) {
        return Err(eyre::eyre!("Path {} not found in graph", node.path_id()));
      }
      let path = &self.paths[&node.path_id()];

      if !block.alignments().contains_key(node_id) {
        return Err(eyre::eyre!("Node {} not found in block {}", node_id, block.id()));
      }

      if !path.nodes.contains(node_id) {
        return Err(eyre::eyre!("Node {} not found in path {}", node_id, path.id()));
      }
    }

    for (block_id, block) in &self.blocks {
      if block.alignments().is_empty() {
        return Err(eyre::eyre!("Block {} has no nodes", block_id));
      }

      for node_id in block.alignments().keys() {
        if !self.nodes.contains_key(node_id) {
          return Err(eyre::eyre!("Node {} not found in graph", node_id));
        }
      }
    }

    for (path_id, path) in &self.paths {
      for node_id in &path.nodes {
        if !self.nodes.contains_key(node_id) {
          return Err(eyre::eyre!("Node {node_id} from path {path_id} not found in graph"));
        }
      }

      // // check that there are no duplicated node ids
      // // currently disabled because this could rarely happen for empty nodes
      // let mut seen = BTreeSet::new();
      // for node_id in &path.nodes {
      //   if !seen.insert(node_id) {
      //     return Err(eyre::eyre!("Node {node_id} appears more than once in path {path_id}",));
      //   }
      // }

      // check that nodes in the same path have contiguous positions
      let mut prev_pos = self.nodes[path.nodes.first().unwrap()].position().1;
      for &node_id in &path.nodes[1..] {
        let pos = self.nodes[&node_id].position().0;
        if pos != prev_pos {
          return Err(eyre::eyre!(
            "Node {node_id} in path {path_id} has position {pos} but previous node has position {prev_pos}",
          ));
        }
        prev_pos = self.nodes[&node_id].position().1;
      }
      if path.circular() {
        let first_pos = self.nodes[path.nodes.first().unwrap()].position().0;
        let last_pos = self.nodes[path.nodes.last().unwrap()].position().1;
        if first_pos != last_pos {
          return Err(eyre::eyre!(
            "Circular path {path_id} has first node position {first_pos} different from last node position {last_pos}",
          ));
        }
      }
    }

    Ok(())
  }

  pub fn path_ids(&self) -> impl Iterator<Item = PathId> + '_ {
    self.paths.keys().copied()
  }

  pub fn block_ids(&self) -> impl Iterator<Item = BlockId> + '_ {
    self.blocks.keys().copied()
  }

  pub fn node_ids(&self) -> impl Iterator<Item = NodeId> + '_ {
    self.nodes.keys().copied()
  }

  pub fn paths(&self) -> impl Iterator<Item = &PangraphPath> {
    self.paths.values()
  }

  pub fn path_names(&self) -> impl Iterator<Item = Option<&str>> {
    self.paths.values().map(|path| path.name.as_deref())
  }

  /// Returns a list of core block ids. Core blocks are present exactly once in each path.
  pub fn core_block_ids(&self) -> impl Iterator<Item = BlockId> + '_ {
    let path_ids: BTreeSet<_> = self.path_ids().collect();
    self.blocks.iter().filter_map(move |(block_id, block)| {
      let block_path_ids: BTreeSet<_> = block
        .alignment_keys()
        .into_iter()
        .map(|nid| self.nodes[&nid].path_id())
        .collect();

      // n. of nodes in the block
      let n_nodes = block.alignment_keys().len();

      // check that the block is present in all paths
      let is_in_all_paths = block_path_ids == path_ids;
      // check that the block is not duplicated in any path,
      // i.e. the number of nodes is equal to the number of path it is present in
      let is_not_duplicated = n_nodes == block_path_ids.len();
      (is_in_all_paths && is_not_duplicated).then_some(*block_id)
    })
  }

  // Returns the path id given the path name.
  pub fn path_id_by_name(&self, path_name: impl AsRef<str>) -> Result<PathId, Report> {
    let path_name = path_name.as_ref();
    self
      .paths
      .iter()
      .find(|(_, path)| path.name.as_deref() == Some(path_name))
      .map(|(pid, _)| *pid)
      .ok_or_else(|| make_internal_report!("When retrieving path id by name: path '{path_name}' not found"))
  }
}

impl FromStr for Pangraph {
  type Err = Report;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    json_read_str(s).wrap_err("When parsing Pangraph JSON contents")
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
      BlockId(1) => PangraphBlock::new(BlockId(1), "1",
        btreemap!{ NodeId(1) => Edit::empty(), NodeId(2) => Edit::empty() }),
      BlockId(2) => PangraphBlock::new(BlockId(2), "2",
        btreemap!{ NodeId(3) => Edit::empty(), NodeId(4) => Edit::empty(), NodeId(5) => Edit::empty() }),
      BlockId(3) => PangraphBlock::new(BlockId(3), "3",
        btreemap!{ NodeId(6) => Edit::empty(), NodeId(7) => Edit::empty(), NodeId(8) => Edit::empty(), }),
    };

    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1), NodeId(3), NodeId(6)], 0, false, None, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [NodeId(4), NodeId(7)           ], 0, false, None, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [NodeId(2), NodeId(5), NodeId(8)], 0, false, None, None),
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
      BlockId(4) => PangraphBlock::new(BlockId(4), "4", btreemap!{}),
      BlockId(5) => PangraphBlock::new(BlockId(5), "5", btreemap!{}),
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
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1),  NodeId(9),  NodeId(10),  NodeId(6)], 0, false, None, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [NodeId(11), NodeId(12), NodeId(7)             ], 0, false, None, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [NodeId(2),  NodeId(14), NodeId(13),  NodeId(8)], 0, false, None, None),
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
