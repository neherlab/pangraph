use crate::io::fasta::FastaRecord;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_read_str;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::pangraph::pangraph_path::{PangraphPath, PathId, genome_seed};
use crate::pangraph::strand::Strand;
use crate::representation::seq::Seq;
use crate::tree::clade::WithNewickName;
use crate::utils::map_merge::{ConflictResolution, map_merge};
use crate::{make_error, make_internal_report, make_report};
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
    // Block and node ids are seeded from the genome name rather than from the record index, so that
    // two graphs built independently cannot collide: names are unique within a build, and `merge`
    // rejects graphs that share one. Seeding from the index instead made every build start at
    // `0, 1, 2, ...`, so merging two graphs required relabeling one of them first.
    //
    // Path ids stay sequential: they double as the ordering index of the genomes, and keying the
    // `paths` map by a hash would scramble genome order in every output.
    let seed = genome_seed(&fasta.seq_name);
    let node_id = NodeId(seed);
    let block_id = BlockId(seed);
    let block = PangraphBlock::from_consensus(fasta.seq, block_id, node_id);
    let path_id = PathId(fasta.index);
    let node_position = if circular { (0, 0) } else { (0, tot_len) }; // path wraps around if circular
    let node = PangraphNode::new(node_id, block.id(), path_id, strand, node_position);
    let path = PangraphPath::new(
      path_id,
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

  /// Reads a graph from a JSON file, or from stdin when no path is given, and validates it.
  ///
  /// The single entry point through which every command reads a graph, and therefore the one place
  /// that has both the graph and the name of the file it came from. [`Self::validate`] runs here so
  /// that a malformed file is rejected up front, named, and before any work is done on it, rather
  /// than surfacing much later as a failed lookup deep inside an algorithm.
  pub fn from_path<P: AsRef<Path>>(filepath: &Option<P>) -> Result<Self, Report> {
    let reader = open_file_or_stdin(filepath)?;
    let data = read_reader_to_string(reader).wrap_err("When reading Pangraph JSON")?;
    let graph = Self::from_str(&data).wrap_err("When parsing Pangraph JSON")?;
    graph.validate().wrap_err_with(|| match filepath {
      Some(filepath) => format!("When validating the graph read from '{}'", filepath.as_ref().display()),
      None => "When validating the graph read from standard input".to_owned(),
    })?;
    Ok(graph)
  }

  pub fn to_string_pretty(&self) -> Result<String, Report> {
    let mut tree_str = serde_json::to_string_pretty(self)?;
    tree_str += "\n";
    Ok(tree_str)
  }

  pub fn consensuses(&self) -> impl Iterator<Item = &Seq> {
    self.blocks.values().map(|block| block.consensus())
  }

  /// Returns this graph with its path ids renumbered contiguously from `offset`, preserving their
  /// relative order.
  ///
  /// Blocks are untouched and node ids are preserved: only the `path_id` a node stores, and the id
  /// of the path itself, are rewritten. Node ids are derived from [`PangraphPath::seed`] rather
  /// than from the path id (see [`PangraphNode::with_derived_id`]), so renumbering cannot put a
  /// node id out of step with the contents it was derived from.
  ///
  /// Used by `merge` to lift one graph's genomes above the other's: path ids are sequential within
  /// each graph, so two graphs built independently always collide on them. The graph is taken by
  /// value so that the sequences can be moved rather than copied.
  ///
  /// Errors if a node refers to a path the graph does not contain, which can only happen in a graph
  /// that pangraph did not write.
  pub fn renumber_paths(self, offset: usize) -> Result<Self, Report> {
    let Self {
      paths,
      blocks,
      mut nodes,
    } = self;

    let path_map: BTreeMap<PathId, PathId> = paths
      .keys()
      .enumerate()
      .map(|(rank, &pid)| (pid, PathId(offset + rank)))
      .collect();

    // Node ids do not change, so the map is updated in place rather than rebuilt: only the path a
    // node points at moves.
    //
    // `path_map` is keyed by the map keys of this graph, so looking a path's own new id up by its
    // key cannot miss. The path id read off a node's *contents* is a different matter: it can
    // dangle in a graph that pangraph did not write, so that lookup is fallible.
    for (nid, node) in &mut nodes {
      let path_id = *path_map.get(&node.path_id()).ok_or_else(|| {
        make_report!(
          "Node {nid} refers to path {}, which the graph does not contain",
          node.path_id()
        )
      })?;
      node.set_path_id(path_id);
    }

    let paths: BTreeMap<PathId, PangraphPath> = paths
      .into_iter()
      .map(|(pid, mut path)| {
        path.id = path_map[&pid];
        (path.id, path)
      })
      .collect();

    Ok(Self { paths, blocks, nodes })
  }

  /// Returns true if this graph shares no block, node or path id with `other`.
  ///
  /// Block and node ids are derived from genome names, which `merge` requires to be distinct across
  /// the two graphs, so this holds by construction for graphs written by pangraph 1.4 or later.
  /// Graphs written by earlier versions derive their ids from the input order instead, and two of
  /// those do collide, which is what `merge` uses this to detect.
  pub fn is_id_disjoint_from(&self, other: &Self) -> bool {
    self.blocks.keys().all(|bid| !other.blocks.contains_key(bid))
      && self.nodes.keys().all(|nid| !other.nodes.contains_key(nid))
      && self.paths.keys().all(|pid| !other.paths.contains_key(pid))
  }

  /// Returns the smallest path id that is above every path id of this graph.
  ///
  /// The offset `merge` renumbers the appended graph from, so that the two sets of path ids cannot
  /// overlap. Not simply the number of paths: `simplify` drops paths without renumbering, so path
  /// ids are not necessarily contiguous.
  pub fn path_id_upper_bound(&self) -> usize {
    self.paths.keys().map(|pid| pid.0 + 1).max().unwrap_or(0)
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
      ConflictResolution::Custom(|(kl, _vl), (_kr, _vr)| panic!("Conflicting key: '{kl}'")),
    );

    for (old_node_id, new_nodes) in &u.n_new {
      let path_id = self.nodes[old_node_id].path_id();

      let path = self.paths.get_mut(&path_id).unwrap(); // FIXME

      let old_idx = path.nodes.iter().position(|node_id| node_id == old_node_id).unwrap();

      // only one such nodes
      debug_assert_eq!(path.nodes.iter().filter(|node_id| *node_id == old_node_id).count(), 1);

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
      .filter(|(_, block)| block.alignments().is_empty())
      .map(|(bid, _)| *bid)
      .collect();

    for bid in empty_blocks {
      self.blocks.remove(&bid);
    }
  }

  /// Checks the invariants that let the rest of pangraph resolve graph entities by id without
  /// error handling, and that keep sequence reconstruction inside array bounds.
  ///
  /// Unlike [`Self::sanity_check`] this runs in release builds, because it guards against a
  /// malformed *input*: a graph read from a file was not necessarily written by pangraph, so
  /// nothing guarantees that its cross-references resolve or that its offsets are in range. Every
  /// graph read through [`Self::from_path`] is validated, and that is what lets the code downstream
  /// index the maps directly and treat a lookup that fails anyway as an internal error.
  ///
  /// Deliberately limited to what a bad file can break. It does *not* check that the graph is
  /// semantically coherent (that node positions tile the genome, that edits do not overlap), since
  /// those are symptoms of a bug in pangraph rather than of a bad input, and are covered by
  /// [`Self::sanity_check`] in debug builds. Errors are reported as ordinary user-facing errors for
  /// the same reason: the offending graph came from the user.
  ///
  /// Costs roughly a tenth of the JSON parse it follows, so it is always worth running.
  pub fn validate(&self) -> Result<(), Report> {
    // Each entity is stored under a map key *and* carries its own id. Serde keys the maps by the
    // JSON object key, so the two can disagree in a graph that pangraph did not write. Everything
    // that resolves an entity by key while reading its id off the entity itself depends on them
    // agreeing, so check it first. One integer comparison per entity.
    for (block_id, block) in &self.blocks {
      if block.id() != *block_id {
        return make_error!("Block is stored under id {block_id} but reports id {}", block.id());
      }
    }

    for (node_id, node) in &self.nodes {
      if node.id() != *node_id {
        return make_error!("Node is stored under id {node_id} but reports id {}", node.id());
      }
    }

    for (path_id, path) in &self.paths {
      if path.id() != *path_id {
        return make_error!("Path is stored under id {path_id} but reports id {}", path.id());
      }
    }

    // Which path walks each node, as claimed by the paths. Collected in one pass because
    // `path.nodes` is a `Vec`: scanning it per node would re-read a genome's entire walk once for
    // every node of that genome.
    let mut walked_by: BTreeMap<NodeId, PathId> = BTreeMap::new();
    for (path_id, path) in &self.paths {
      for node_id in &path.nodes {
        if !self.nodes.contains_key(node_id) {
          return make_error!("Node {node_id} from path {path_id} not found in graph");
        }
        // A node id may repeat within the walk of one path, which happens for empty nodes, but two
        // paths sharing a node would make the node's own `path_id` ambiguous.
        if let Some(previous) = walked_by.insert(*node_id, *path_id) {
          if previous != *path_id {
            return make_error!("Node {node_id} is walked by both path {previous} and path {path_id}");
          }
        }
      }
    }

    for (node_id, node) in &self.nodes {
      let Some(block) = self.blocks.get(&node.block_id()) else {
        return make_error!("Block {} of node {node_id} not found in graph", node.block_id());
      };
      if !block.alignments().contains_key(node_id) {
        return make_error!("Node {node_id} not found in block {}", block.id());
      }

      let Some(path) = self.paths.get(&node.path_id()) else {
        return make_error!("Path {} of node {node_id} not found in graph", node.path_id());
      };
      if walked_by.get(node_id) != Some(&node.path_id()) {
        return make_error!("Node {node_id} is not in the walk of its path {}", node.path_id());
      }

      // Reconstruction rotates a genome by the offset of its first node, so an offset reaching past
      // the end of that genome would rotate by more than the genome's length. Both ends are
      // compared, and `tot_len` itself is a valid offset: the last node of a circular path ends
      // where the genome does.
      let (start, end) = node.position();
      if start > path.tot_len() || end > path.tot_len() {
        return make_error!(
          "Node {node_id} has position ({start}, {end}), outside its path {} of total length {}",
          path.id(),
          path.tot_len()
        );
      }
    }

    for (block_id, block) in &self.blocks {
      let consensus_len = block.consensus().len();
      for (node_id, edits) in block.alignments() {
        if !self.nodes.contains_key(node_id) {
          return make_error!("Node {node_id} of block {block_id} not found in graph");
        }
        // `Edit::apply` indexes the consensus by edit position while reconstructing this node.
        edits
          .check_bounds(consensus_len)
          .wrap_err_with(|| format!("When checking the alignment of node {node_id} against block {block_id}"))?;
      }
    }

    Ok(())
  }

  #[cfg(any(test, debug_assertions))]
  pub fn sanity_check(&self) -> Result<(), Report> {
    // Referential integrity and the bounds that keep reconstruction in range. Shared with the
    // release-build validation of graphs read from a file, since a graph that pangraph just built
    // must satisfy at least as much as one it is willing to load.
    self.validate()?;

    for (block_id, block) in &self.blocks {
      if block.alignments().is_empty() {
        return Err(eyre::eyre!("Block {} has no nodes", block_id));
      }
    }

    for (path_id, path) in &self.paths {
      // // check that there are no duplicated node ids
      // // currently disabled because this could rarely happen for empty nodes
      // let mut seen = BTreeSet::new();
      // for node_id in &path.nodes {
      //   if !seen.insert(node_id) {
      //     return Err(eyre::eyre!("Node {node_id} appears more than once in path {path_id}",));
      //   }
      // }

      // check that nodes in the same path have contiguous positions
      if let Some(first_node_id) = path.nodes.first() {
        let mut prev_pos = self.nodes[first_node_id].position().1;
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
          let first_pos = self.nodes[first_node_id].position().0;
          let last_node_id = path
            .nodes
            .last()
            .ok_or_else(|| eyre::eyre!("Path {path_id} has first node but no last node"))?;
          let last_pos = self.nodes[last_node_id].position().1;
          if first_pos != last_pos {
            return Err(eyre::eyre!(
              "Circular path {path_id} has first node position {first_pos} different from last node position {last_pos}",
            ));
          }
        }
      }
    }

    Ok(())
  }

  pub fn path_ids(&self) -> impl Iterator<Item = PathId> + use<'_> {
    self.paths.keys().copied()
  }

  pub fn block_ids(&self) -> impl Iterator<Item = BlockId> + use<'_> {
    self.blocks.keys().copied()
  }

  pub fn node_ids(&self) -> impl Iterator<Item = NodeId> + use<'_> {
    self.nodes.keys().copied()
  }

  pub fn paths(&self) -> impl Iterator<Item = &PangraphPath> {
    self.paths.values()
  }

  pub fn path_names(&self) -> impl Iterator<Item = Option<&str>> {
    self.paths.values().map(|path| path.name.as_deref())
  }

  /// Returns a list of core block ids. Core blocks are present exactly once in each path.
  pub fn core_block_ids(&self) -> impl Iterator<Item = BlockId> + use<'_> {
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

impl WithNewickName for Pangraph {
  /// Returns a Newick-safe label for this graph: the path name for a singleton,
  /// or names joined with `|` for a multi-path graph. `None` if no paths are named.
  fn newick_name(&self) -> Option<String> {
    let names: Vec<&str> = self.paths.values().filter_map(|p| p.name.as_deref()).collect();
    (!names.is_empty()).then(|| names.join("|"))
  }
}

impl WithNewickName for Option<Pangraph> {
  /// Internal tree nodes carry `None` (unlabeled in Newick); leaves delegate to `Pangraph`.
  fn newick_name(&self) -> Option<String> {
    self.as_ref().and_then(WithNewickName::newick_name)
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
  use crate::o;
  use crate::pangraph::edits::{Edit, Sub};
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::pangraph_path::PangraphPath;
  use crate::pangraph::reconstruct::reconstruct;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::utils::error::report_to_string;
  use itertools::Itertools;
  use maplit::btreemap;
  use rstest::rstest;

  #[test]
  fn test_graph_update() {
    // graph
    // p1 -> [b1+,b2+,b3+]
    // p2 -> [b2+,b3+]
    // p3 -> [b1+,b2-,b3+]
    // update
    // b2+ -> [b4+, b5-]

    let nodes = btreemap! {
      NodeId(1) => PangraphNode::new(NodeId(1), BlockId(1), PathId(1), Forward,  (0, 0)), // FIXME
      NodeId(2) => PangraphNode::new(NodeId(2), BlockId(1), PathId(3), Forward,  (0, 0)), // FIXME
      NodeId(3) => PangraphNode::new(NodeId(3), BlockId(2), PathId(1), Forward,  (0, 0)), // FIXME
      NodeId(4) => PangraphNode::new(NodeId(4), BlockId(2), PathId(2), Forward,  (0, 0)), // FIXME
      NodeId(5) => PangraphNode::new(NodeId(5), BlockId(2), PathId(3), Reverse, (0, 0)), // FIXME
      NodeId(6) => PangraphNode::new(NodeId(6), BlockId(3), PathId(1), Forward,  (0, 0)), // FIXME
      NodeId(7) => PangraphNode::new(NodeId(7), BlockId(3), PathId(2), Forward,  (0, 0)), // FIXME
      NodeId(8) => PangraphNode::new(NodeId(8), BlockId(3), PathId(3), Forward,  (0, 0)) // FIXME
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
      PathId(1) => PangraphPath::new(PathId(1), [NodeId(1), NodeId(3), NodeId(6)], 0, false, None, None),
      PathId(2) => PangraphPath::new(PathId(2), [NodeId(4), NodeId(7)           ], 0, false, None, None),
      PathId(3) => PangraphPath::new(PathId(3), [NodeId(2), NodeId(5), NodeId(8)], 0, false, None, None),
    };

    let mut G = Pangraph {
      paths: paths.clone(),
      blocks: blocks.clone(),
      nodes: nodes.clone(),
    };

    let new_nodes = btreemap! {
      NodeId(9)  => PangraphNode::new(NodeId(9),  BlockId(4), PathId(1), Forward,  (0, 0)),
      NodeId(10) => PangraphNode::new(NodeId(10), BlockId(5), PathId(1), Reverse, (0, 0)),
      NodeId(11) => PangraphNode::new(NodeId(11), BlockId(4), PathId(2), Forward,  (0, 0)),
      NodeId(12) => PangraphNode::new(NodeId(12), BlockId(5), PathId(2), Reverse, (0, 0)),
      NodeId(13) => PangraphNode::new(NodeId(13), BlockId(4), PathId(3), Reverse, (0, 0)),
      NodeId(14) => PangraphNode::new(NodeId(14), BlockId(5), PathId(3), Forward,  (0, 0)),
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
      PathId(1) => PangraphPath::new(PathId(1), [NodeId(1),  NodeId(9),  NodeId(10),  NodeId(6)], 0, false, None, None),
      PathId(2) => PangraphPath::new(PathId(2), [NodeId(11), NodeId(12), NodeId(7)             ], 0, false, None, None),
      PathId(3) => PangraphPath::new(PathId(3), [NodeId(2),  NodeId(14), NodeId(13),  NodeId(8)], 0, false, None, None),
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

  /// Builds a `Pangraph` whose only relevant content for `newick_name` is its `paths` map.
  /// Blocks and nodes are left empty since `newick_name` only inspects path names.
  fn pangraph_with_named_paths(names: &[Option<&str>]) -> Pangraph {
    let paths = names
      .iter()
      .enumerate()
      .map(|(i, name)| {
        let path = PangraphPath::new(PathId(i), Vec::<NodeId>::new(), 0, false, name.map(String::from), None);
        (path.id, path)
      })
      .collect::<BTreeMap<_, _>>();
    Pangraph {
      paths,
      blocks: BTreeMap::new(),
      nodes: BTreeMap::new(),
    }
  }

  #[test]
  fn test_newick_no_graph() {
    let g: Option<Pangraph> = None;
    assert_eq!(g.newick_name(), None);
  }

  #[rstest]
  #[case::singleton_named(&[Some("isolate_A")], Some("isolate_A".to_owned()))]
  #[case::singleton_unnamed(&[None], None)]
  #[case::multi_path_all_named(&[Some("a"), Some("b"), Some("c")], Some("a|b|c".to_owned()))]
  #[case::multi_path_some_unnamed(&[Some("a"), None, Some("c")], Some("a|c".to_owned()))]
  fn test_newick_name(#[case] names: &[Option<&str>], #[case] expected: Option<String>) {
    let g = pangraph_with_named_paths(names);
    assert_eq!(g.newick_name(), expected);
  }

  /// Builds a two-genome graph the way `build` would, with block and node ids seeded from the
  /// genome names. Two such graphs collide only on their path ids.
  fn two_genome_graph(names: [&str; 2]) -> Pangraph {
    let seeds = names.map(genome_seed);
    let (b0, b1) = (BlockId(seeds[0]), BlockId(seeds[1]));
    let (n0, n1) = (NodeId(seeds[0]), NodeId(seeds[1]));

    let blocks = btreemap! {
      b0 => PangraphBlock::new(b0, "ACGTACGT", btreemap!{ n0 => Edit::empty() }),
      b1 => PangraphBlock::new(b1, "TTTTGGGG", btreemap!{ n1 => Edit::empty() }),
    };
    let nodes = btreemap! {
      n0 => PangraphNode::new(n0, b0, PathId(0), Forward, (0, 8)),
      n1 => PangraphNode::new(n1, b1, PathId(1), Reverse, (0, 8)),
    };
    let paths = btreemap! {
      PathId(0) => PangraphPath::new(PathId(0), [n0], 8, false, Some(names[0].to_owned()), None),
      PathId(1) => PangraphPath::new(PathId(1), [n1], 8, false, Some(names[1].to_owned()), None),
    };
    Pangraph { paths, blocks, nodes }
  }

  /// A graph read from JSON is keyed by the object key, while each entity also stores its own id;
  /// serde never checks that the two agree. `renumber_paths` resolves paths by key, so a
  /// disagreement would otherwise abort the process with a bare "no entry found for key" panic and
  /// no indication of which file was at fault.
  #[rstest]
  fn test_renumber_paths_reports_node_referring_to_a_missing_path() {
    let mut graph = two_genome_graph(["a", "b"]);
    let nid = graph.node_ids().next().unwrap();
    let node = &graph.nodes[&nid];
    let detached = PangraphNode::new(nid, node.block_id(), PathId(99), node.strand(), node.position());
    graph.nodes.insert(nid, detached);

    let err = report_to_string(&graph.renumber_paths(7).unwrap_err());
    assert!(
      err.contains("refers to path 99, which the graph does not contain"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_renumber_paths_keeps_graph_consistent() {
    let original = two_genome_graph(["a", "b"]);
    let graph = original.clone().renumber_paths(7).unwrap();

    graph.sanity_check().unwrap();
    assert_eq!(graph.blocks.len(), 2);
    assert_eq!(graph.nodes.len(), 2);
    assert_eq!(graph.paths.len(), 2);

    // path ids are renumbered contiguously from the offset, in their original order
    assert_eq!(graph.path_ids().collect_vec(), vec![PathId(7), PathId(8)]);
    assert_eq!(
      graph.paths.values().map(|p| p.name.clone()).collect_vec(),
      vec![Some(o!("a")), Some(o!("b"))]
    );

    // blocks and node ids are untouched: only the path id a node stores is rewritten
    assert_eq!(graph.blocks, original.blocks);
    assert_eq!(graph.node_ids().collect_vec(), original.node_ids().collect_vec());
    for path in graph.paths.values() {
      for nid in path.nodes() {
        assert_eq!(graph.nodes[nid].path_id(), path.id());
      }
    }
  }

  #[rstest]
  fn test_renumber_paths_preserves_sequences() {
    let original = two_genome_graph(["a", "b"]);
    let renumbered = original.clone().renumber_paths(5).unwrap();

    let seqs = |g: &Pangraph| {
      reconstruct(g)
        .map(|r| r.map(|r| (r.seq_name, r.seq)))
        .collect::<Result<BTreeMap<_, _>, Report>>()
        .unwrap()
    };

    assert_eq!(seqs(&original), seqs(&renumbered));
  }

  /// `simplify` drops paths without renumbering the survivors, so the offset cannot be taken to be
  /// the number of paths.
  #[rstest]
  fn test_path_id_upper_bound_clears_sparse_path_ids() {
    let mut graph = two_genome_graph(["a", "b"]);
    assert_eq!(graph.path_id_upper_bound(), 2);

    graph.paths.remove(&PathId(0));
    assert_eq!(graph.paths.len(), 1);
    assert_eq!(graph.path_id_upper_bound(), 2);
  }

  #[rstest]
  fn test_path_id_upper_bound_of_an_empty_graph() {
    let graph = Pangraph {
      paths: btreemap! {},
      blocks: btreemap! {},
      nodes: btreemap! {},
    };
    assert_eq!(graph.path_id_upper_bound(), 0);
  }

  /// The property the whole identifier model rests on: two graphs built from differently named
  /// genomes share no block or node id, so `merge` can join them without relabeling either.
  #[rstest]
  fn test_graphs_with_distinct_genome_names_are_id_disjoint() {
    let left = two_genome_graph(["a", "b"]);
    let right = two_genome_graph(["c", "d"])
      .renumber_paths(left.path_id_upper_bound())
      .unwrap();

    assert!(right.is_id_disjoint_from(&left));
    right.sanity_check().unwrap();

    assert_eq!(left.path_ids().collect_vec(), vec![PathId(0), PathId(1)]);
    assert_eq!(right.path_ids().collect_vec(), vec![PathId(2), PathId(3)]);

    // joining the two does not conflict
    let joined = crate::pangraph::graph_merging::graph_join(&left, &right);
    joined.sanity_check().unwrap();
    assert_eq!(joined.paths.len(), 4);
    assert_eq!(joined.blocks.len(), 4);
    assert_eq!(joined.nodes.len(), 4);
  }

  /// Appending to a graph that already absorbed another one. This needed a varying relabeling salt
  /// back when block and node ids were re-derived on merge; deriving them from genome names instead
  /// makes it hold with no bookkeeping at all.
  #[rstest]
  fn test_appending_to_an_already_merged_graph() {
    let first = two_genome_graph(["a", "b"]);
    let second = two_genome_graph(["c", "d"])
      .renumber_paths(first.path_id_upper_bound())
      .unwrap();
    let joined = crate::pangraph::graph_merging::graph_join(&first, &second);

    let third = two_genome_graph(["e", "f"])
      .renumber_paths(joined.path_id_upper_bound())
      .unwrap();

    assert!(third.is_id_disjoint_from(&joined));
    third.sanity_check().unwrap();
    assert_eq!(third.path_ids().collect_vec(), vec![PathId(4), PathId(5)]);

    // and the three of them can be joined without conflicts
    let joined = crate::pangraph::graph_merging::graph_join(&joined, &third);
    joined.sanity_check().unwrap();
    assert_eq!(joined.paths.len(), 6);
    assert_eq!(joined.blocks.len(), 6);
    assert_eq!(joined.nodes.len(), 6);
  }

  /// Two graphs holding the same genome do collide, which is why `merge` rejects them by name
  /// before it ever gets as far as joining them.
  #[rstest]
  fn test_graphs_sharing_a_genome_name_are_not_id_disjoint() {
    let left = two_genome_graph(["a", "b"]);
    let right = two_genome_graph(["a", "c"])
      .renumber_paths(left.path_id_upper_bound())
      .unwrap();

    assert!(!right.is_id_disjoint_from(&left));
  }

  /// Ids no longer depend on the order the input sequences were read in, only on their names.
  #[rstest]
  fn test_singleton_ids_are_seeded_from_the_name_not_the_index() {
    let singleton = |name: &str, index: usize| {
      Pangraph::singleton(
        FastaRecord {
          seq_name: name.to_owned(),
          desc: None,
          seq: Seq::from_str("ACGTACGT"),
          index,
        },
        Forward,
        false,
      )
    };

    // the same genome read at a different position in the input gets the same block and node ids
    let first = singleton("a", 0);
    let shifted = singleton("a", 7);
    assert_eq!(first.block_ids().collect_vec(), shifted.block_ids().collect_vec());
    assert_eq!(first.node_ids().collect_vec(), shifted.node_ids().collect_vec());

    // but the path id still records where it was read, so genome order survives
    assert_eq!(first.path_ids().collect_vec(), vec![PathId(0)]);
    assert_eq!(shifted.path_ids().collect_vec(), vec![PathId(7)]);

    // and a different genome read at the same position gets different block and node ids
    let other = singleton("b", 0);
    assert_ne!(first.block_ids().collect_vec(), other.block_ids().collect_vec());
    assert_ne!(first.node_ids().collect_vec(), other.node_ids().collect_vec());
  }

  /// `validate` is the contract every graph read from a file must satisfy, and the reason the code
  /// downstream is allowed to resolve ids by direct indexing. Each case below is a malformation
  /// that a hand-edited or third-party JSON graph can carry, and that used to reach an indexing
  /// panic in a release build, where `sanity_check` is compiled out.
  #[rstest]
  fn test_validate_accepts_a_well_formed_graph() {
    two_genome_graph(["a", "b"]).validate().unwrap();
  }

  #[rstest]
  fn test_validate_rejects_path_referring_to_a_missing_node() {
    let mut graph = two_genome_graph(["a", "b"]);
    graph.paths.get_mut(&PathId(0)).unwrap().nodes = vec![NodeId(99)];

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Node 99 from path 0 not found in graph"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_validate_rejects_node_referring_to_a_missing_block() {
    let mut graph = two_genome_graph(["a", "b"]);
    let nid = graph.paths[&PathId(0)].nodes[0];
    let node = &graph.nodes[&nid];
    graph.nodes.insert(
      nid,
      PangraphNode::new(nid, BlockId(99), node.path_id(), node.strand(), node.position()),
    );

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains(&format!("Block 99 of node {nid} not found in graph")),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_validate_rejects_node_referring_to_a_missing_path() {
    let mut graph = two_genome_graph(["a", "b"]);
    let nid = graph.paths[&PathId(0)].nodes[0];
    let node = &graph.nodes[&nid];
    graph.nodes.insert(
      nid,
      PangraphNode::new(nid, node.block_id(), PathId(99), node.strand(), node.position()),
    );

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains(&format!("Path 99 of node {nid} not found in graph")),
      "unexpected error: {err}"
    );
  }

  /// A node whose path exists but does not walk it: reconstruction would never emit this node, and
  /// its `path_id` claims a genome it is not part of.
  #[rstest]
  fn test_validate_rejects_node_missing_from_the_walk_of_its_path() {
    let mut graph = two_genome_graph(["a", "b"]);
    let nid = graph.paths[&PathId(0)].nodes[0];
    graph.paths.get_mut(&PathId(0)).unwrap().nodes = vec![];

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains(&format!("Node {nid} is not in the walk of its path 0")),
      "unexpected error: {err}"
    );
  }

  /// Node ids are unique across the graph, so a node claimed by two walks makes its own `path_id`
  /// meaningless, and would have the node reconstructed into two different genomes.
  #[rstest]
  fn test_validate_rejects_node_walked_by_two_paths() {
    let mut graph = two_genome_graph(["a", "b"]);
    let (first, second) = (graph.paths[&PathId(0)].nodes[0], graph.paths[&PathId(1)].nodes[0]);
    graph.paths.get_mut(&PathId(1)).unwrap().nodes = vec![first, second];

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains(&format!("Node {first} is walked by both path 0 and path 1")),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_validate_rejects_node_missing_from_the_alignments_of_its_block() {
    let mut graph = two_genome_graph(["a", "b"]);
    let nid = graph.paths[&PathId(0)].nodes[0];
    let bid = graph.nodes[&nid].block_id();
    graph
      .blocks
      .insert(bid, PangraphBlock::new(bid, "ACGTACGT", btreemap! {}));

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains(&format!("Node {nid} not found in block {bid}")),
      "unexpected error: {err}"
    );
  }

  /// Reconstruction rotates a genome by the offset of its first node. An offset past the end of the
  /// genome used to reach `rotate_right`, which panics rather than reporting.
  #[rstest]
  fn test_validate_rejects_node_position_beyond_the_length_of_its_path() {
    let mut graph = two_genome_graph(["a", "b"]);
    let nid = graph.paths[&PathId(0)].nodes[0];
    let node = &graph.nodes[&nid];
    graph.nodes.insert(
      nid,
      PangraphNode::new(nid, node.block_id(), node.path_id(), node.strand(), (100, 8)),
    );

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains(&format!(
        "Node {nid} has position (100, 8), outside its path 0 of total length 8"
      )),
      "unexpected error: {err}"
    );
  }

  /// `Edit::apply` indexes the consensus by edit position, so an out-of-range edit used to panic
  /// while reconstructing the block, even though `apply` returns a `Result`.
  #[rstest]
  fn test_validate_rejects_edit_beyond_the_consensus_of_its_block() {
    let mut graph = two_genome_graph(["a", "b"]);
    let nid = graph.paths[&PathId(0)].nodes[0];
    let bid = graph.nodes[&nid].block_id();
    let edit = Edit {
      subs: vec![Sub::new(99, 'T')],
      dels: vec![],
      inss: vec![],
    };
    graph
      .blocks
      .insert(bid, PangraphBlock::new(bid, "ACGTACGT", btreemap! { nid => edit }));

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Substitution position 99 is out of bounds for sequence of length 8"),
      "unexpected error: {err}"
    );
    assert!(
      err.contains(&format!("alignment of node {nid} against block {bid}")),
      "unexpected error: {err}"
    );
  }
}
