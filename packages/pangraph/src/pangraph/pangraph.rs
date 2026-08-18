use crate::io::fasta::FastaRecord;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_read_str;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::pangraph::strand::Strand;
use crate::representation::seq::Seq;
use crate::tree::clade::WithNewickName;
use crate::utils::id::id;
use crate::utils::map_merge::{ConflictResolution, map_merge};
use crate::{make_error, make_internal_error, make_internal_report, make_report};
use eyre::{Report, WrapErr};
use log::warn;
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

  /// Returns this graph with every block, node and path id re-derived.
  ///
  /// Block and node ids become `id((salt, old_id))`; path ids are renumbered contiguously from
  /// `path_id_offset`, preserving their relative order (path ids double as the ordering index of
  /// the genomes, so they are kept small and sequential rather than hashed).
  ///
  /// Consensuses, edits, names, descriptions, strands and positions are moved over untouched: the
  /// relabeled graph describes exactly the same sequences as the original one. The graph is taken
  /// by value so that the sequences can be moved rather than copied.
  pub fn relabel(self, salt: usize, path_id_offset: usize) -> Result<Self, Report> {
    let Self { paths, blocks, nodes } = self;
    let (n_blocks, n_nodes, n_paths) = (blocks.len(), nodes.len(), paths.len());

    let block_map: BTreeMap<BlockId, BlockId> = blocks.keys().map(|&bid| (bid, BlockId(id((salt, bid))))).collect();

    let node_map: BTreeMap<NodeId, NodeId> = nodes.keys().map(|&nid| (nid, NodeId(id((salt, nid))))).collect();

    let path_map: BTreeMap<PathId, PathId> = paths
      .keys()
      .enumerate()
      .map(|(rank, &pid)| (pid, PathId(path_id_offset + rank)))
      .collect();

    // The three maps above are keyed by the *map keys* of this graph, so looking an entity's own new
    // id up by its key cannot miss. Ids taken from an entity's *contents* are a different matter:
    // they can dangle, or disagree with the key they are stored under, in a graph that pangraph did
    // not write, so those lookups are fallible.
    let blocks: BTreeMap<BlockId, PangraphBlock> = blocks
      .into_iter()
      .map(|(bid, block)| {
        let block = block.relabel(block_map[&bid], &node_map)?;
        Ok((block.id(), block))
      })
      .collect::<Result<_, Report>>()?;

    let nodes: BTreeMap<NodeId, PangraphNode> = nodes
      .into_iter()
      .map(|(nid, node)| {
        let block_id = *block_map.get(&node.block_id()).ok_or_else(|| {
          make_report!(
            "Node {nid} refers to block {}, which the graph does not contain",
            node.block_id()
          )
        })?;
        let path_id = *path_map.get(&node.path_id()).ok_or_else(|| {
          make_report!(
            "Node {nid} refers to path {}, which the graph does not contain",
            node.path_id()
          )
        })?;

        let node = PangraphNode::new(Some(node_map[&nid]), block_id, path_id, node.strand(), node.position());
        Ok((node.id(), node))
      })
      .collect::<Result<_, Report>>()?;

    let paths: BTreeMap<PathId, PangraphPath> = paths
      .into_iter()
      .map(|(pid, mut path)| {
        path.id = path_map[&pid];
        for nid in &mut path.nodes {
          let old = *nid;
          *nid = *node_map
            .get(&old)
            .ok_or_else(|| make_report!("Path {pid} refers to node {old}, which the graph does not contain"))?;
        }
        Ok((path.id, path))
      })
      .collect::<Result<_, Report>>()?;

    // An injective relabeling cannot change the number of entities. If it did, two distinct ids
    // were mapped onto the same one and entities were silently dropped.
    if (blocks.len(), nodes.len(), paths.len()) != (n_blocks, n_nodes, n_paths) {
      return make_internal_error!(
        "When relabeling graph ids: expected {n_blocks} blocks, {n_nodes} nodes and {n_paths} paths, but got {} blocks, {} nodes and {} paths",
        blocks.len(),
        nodes.len(),
        paths.len(),
      );
    }

    Ok(Self { paths, blocks, nodes })
  }

  /// Returns true if this graph shares no block, node or path id with `other`.
  pub fn is_id_disjoint_from(&self, other: &Self) -> bool {
    self.blocks.keys().all(|bid| !other.blocks.contains_key(bid))
      && self.nodes.keys().all(|nid| !other.nodes.contains_key(nid))
      && self.paths.keys().all(|pid| !other.paths.contains_key(pid))
  }

  /// Relabels this graph so that it shares no identifier with `other`, which is left untouched.
  ///
  /// Two graphs built independently always collide: `Pangraph::singleton` labels the first genome
  /// of every build with path, block and node id `0`, and blocks that never merge keep that id all
  /// the way to the final graph. Merging therefore requires namespacing one of the two graphs
  /// first.
  pub fn make_disjoint_from(self, other: &Self) -> Result<Self, Report> {
    const MAX_ATTEMPTS: usize = 8;

    // Path ids are assigned above every path id of `other`, so they cannot collide by construction.
    let path_id_offset = other.paths.keys().map(|pid| pid.0 + 1).max().unwrap_or(0);

    // The salt must differ from the one used by any earlier merge whose relabeled ids survive into
    // `other`, otherwise those ids get re-derived a second time and land on themselves. The path id
    // offset provides that: a merger salted with `offset(P)` produces a graph with at least one more
    // genome than `P`, so a later merge with that graph on the left uses a strictly larger offset.
    let mut relabeled = self;
    for attempt in 0..MAX_ATTEMPTS {
      relabeled = relabeled.relabel(id((path_id_offset, attempt)), path_id_offset)?;

      if relabeled.is_id_disjoint_from(other) {
        return Ok(relabeled);
      }

      // Relabeling is a composition of injective maps, so retrying on top of the previous attempt is
      // safe, and path ids are assigned by rank rather than derived from the previous id, so they do
      // not drift. Reaching this point needs either a genuine hash collision, or an operation that
      // breaks the monotonicity of the offset: `simplify` drops paths without renumbering, so
      // `build -> merge -> merge -> simplify -> merge` can bring the offset back to a value already
      // used as a salt.
      warn!("Identifier collision when relabeling graph ids (attempt {attempt}); retrying");
    }

    make_internal_error!(
      "When making graphs id-disjoint: no collision-free relabeling of block and node ids found after {MAX_ATTEMPTS} attempts"
    )
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
  /// semantically coherent — that node positions tile the genome, that edits do not overlap — since
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

  /// Builds a `Pangraph` whose only relevant content for `newick_name` is its `paths` map.
  /// Blocks and nodes are left empty since `newick_name` only inspects path names.
  fn pangraph_with_named_paths(names: &[Option<&str>]) -> Pangraph {
    let paths = names
      .iter()
      .enumerate()
      .map(|(i, name)| {
        let path = PangraphPath::new(
          Some(PathId(i)),
          Vec::<NodeId>::new(),
          0,
          false,
          name.map(String::from),
          None,
        );
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

  /// Builds a two-genome graph whose ids are the small sequential integers that `build` assigns to
  /// blocks that never merge. Two such graphs collide on every single id.
  fn colliding_graph(names: [&str; 2]) -> Pangraph {
    let blocks = btreemap! {
      BlockId(0) => PangraphBlock::new(BlockId(0), "ACGTACGT", btreemap!{ NodeId(0) => Edit::empty() }),
      BlockId(1) => PangraphBlock::new(BlockId(1), "TTTTGGGG", btreemap!{ NodeId(1) => Edit::empty() }),
    };
    let nodes = btreemap! {
      NodeId(0) => PangraphNode::new(Some(NodeId(0)), BlockId(0), PathId(0), Forward, (0, 8)),
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(1), PathId(1), Reverse, (0, 8)),
    };
    let paths = btreemap! {
      PathId(0) => PangraphPath::new(Some(PathId(0)), [NodeId(0)], 8, false, Some(names[0].to_owned()), None),
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1)], 8, false, Some(names[1].to_owned()), None),
    };
    Pangraph { paths, blocks, nodes }
  }

  /// A graph read from JSON is keyed by the object key, while each entity also stores its own id;
  /// serde never checks that the two agree. `relabel` resolves entities by key, so a disagreement
  /// used to abort the process with a bare "no entry found for key" panic and no indication of
  /// which file was at fault. It must be a reportable error instead.
  #[rstest]
  fn test_relabel_reports_block_stored_under_a_mismatched_id() {
    let mut graph = colliding_graph(["a", "b"]);
    let block = graph.blocks.remove(&BlockId(0)).unwrap();
    graph.blocks.insert(BlockId(42), block); // key 42, but the block still reports id 0

    let err = report_to_string(&graph.clone().relabel(1, 7).unwrap_err());
    assert!(err.contains("does not contain"), "unexpected error: {err}");

    // `sanity_check` never compared key against stored id, so it used to pass this graph through.
    let err = report_to_string(&graph.sanity_check().unwrap_err());
    assert!(err.contains("stored under id 42"), "unexpected error: {err}");
    assert!(err.contains("reports id 0"), "unexpected error: {err}");
  }

  #[rstest]
  fn test_relabel_reports_path_referring_to_a_missing_node() {
    let mut graph = colliding_graph(["a", "b"]);
    graph.paths.get_mut(&PathId(0)).unwrap().nodes = vec![NodeId(99)];

    let err = report_to_string(&graph.relabel(1, 7).unwrap_err());
    assert!(
      err.contains("Path 0 refers to node 99, which the graph does not contain"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_relabel_reports_node_referring_to_a_missing_block() {
    let mut graph = colliding_graph(["a", "b"]);
    let node = graph.nodes.get_mut(&NodeId(0)).unwrap();
    *node = PangraphNode::new(
      Some(NodeId(0)),
      BlockId(99),
      node.path_id(),
      node.strand(),
      node.position(),
    );

    let err = report_to_string(&graph.relabel(1, 7).unwrap_err());
    assert!(
      err.contains("Node 0 refers to block 99, which the graph does not contain"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_relabel_keeps_graph_consistent() {
    let graph = colliding_graph(["a", "b"]).relabel(1, 7).unwrap();

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

    // block and node ids are re-derived, and no longer the original small integers
    assert!(graph.block_ids().all(|bid| bid.0 > 1));
    assert!(graph.node_ids().all(|nid| nid.0 > 1));
  }

  #[rstest]
  fn test_relabel_preserves_sequences() {
    let original = colliding_graph(["a", "b"]);
    let relabeled = original.clone().relabel(3, 0).unwrap();

    let seqs = |g: &Pangraph| {
      reconstruct(g)
        .map(|r| r.map(|r| (r.seq_name, r.seq)))
        .collect::<Result<BTreeMap<_, _>, Report>>()
        .unwrap()
    };

    assert_eq!(seqs(&original), seqs(&relabeled));
  }

  #[rstest]
  fn test_relabel_is_deterministic() {
    let first = colliding_graph(["a", "b"]).relabel(2, 5).unwrap();
    let second = colliding_graph(["a", "b"]).relabel(2, 5).unwrap();
    assert_eq!(first, second);
  }

  #[rstest]
  fn test_make_disjoint_from() {
    let left = colliding_graph(["a", "b"]);
    let right = colliding_graph(["c", "d"]);

    // the two graphs share every single id before relabeling
    assert!(!right.is_id_disjoint_from(&left));

    let right = right.make_disjoint_from(&left).unwrap();

    assert!(right.is_id_disjoint_from(&left));
    right.sanity_check().unwrap();

    // the left graph is untouched, and the right graph's genomes follow it
    assert_eq!(left.path_ids().collect_vec(), vec![PathId(0), PathId(1)]);
    assert_eq!(right.path_ids().collect_vec(), vec![PathId(2), PathId(3)]);

    // joining the two no longer conflicts
    let joined = crate::pangraph::graph_merging::graph_join(&left, &right);
    joined.sanity_check().unwrap();
    assert_eq!(joined.paths.len(), 4);
    assert_eq!(joined.blocks.len(), 4);
    assert_eq!(joined.nodes.len(), 4);
  }

  /// Appending to a graph that already absorbed a relabeled graph. With a constant salt the ids of
  /// the third graph were re-derived exactly onto those the second one left behind, so the second
  /// append always failed. Every graph here carries the same small ids, which is what `build`
  /// assigns to blocks and nodes that never merge.
  #[rstest]
  fn test_make_disjoint_from_after_a_previous_merge() {
    let first = colliding_graph(["a", "b"]);
    let second = colliding_graph(["c", "d"]).make_disjoint_from(&first).unwrap();
    let joined = crate::pangraph::graph_merging::graph_join(&first, &second);

    let third = colliding_graph(["e", "f"]).make_disjoint_from(&joined).unwrap();

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

  /// `validate` is the contract every graph read from a file must satisfy, and the reason the code
  /// downstream is allowed to resolve ids by direct indexing. Each case below is a malformation
  /// that a hand-edited or third-party JSON graph can carry, and that used to reach an indexing
  /// panic in a release build, where `sanity_check` is compiled out.
  #[rstest]
  fn test_validate_accepts_a_well_formed_graph() {
    colliding_graph(["a", "b"]).validate().unwrap();
  }

  #[rstest]
  fn test_validate_rejects_path_referring_to_a_missing_node() {
    let mut graph = colliding_graph(["a", "b"]);
    graph.paths.get_mut(&PathId(0)).unwrap().nodes = vec![NodeId(99)];

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Node 99 from path 0 not found in graph"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_validate_rejects_node_referring_to_a_missing_block() {
    let mut graph = colliding_graph(["a", "b"]);
    let node = &graph.nodes[&NodeId(0)];
    graph.nodes.insert(
      NodeId(0),
      PangraphNode::new(
        Some(NodeId(0)),
        BlockId(99),
        node.path_id(),
        node.strand(),
        node.position(),
      ),
    );

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Block 99 of node 0 not found in graph"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_validate_rejects_node_referring_to_a_missing_path() {
    let mut graph = colliding_graph(["a", "b"]);
    let node = &graph.nodes[&NodeId(0)];
    graph.nodes.insert(
      NodeId(0),
      PangraphNode::new(
        Some(NodeId(0)),
        node.block_id(),
        PathId(99),
        node.strand(),
        node.position(),
      ),
    );

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Path 99 of node 0 not found in graph"),
      "unexpected error: {err}"
    );
  }

  /// A node whose path exists but does not walk it: reconstruction would never emit this node, and
  /// its `path_id` claims a genome it is not part of.
  #[rstest]
  fn test_validate_rejects_node_missing_from_the_walk_of_its_path() {
    let mut graph = colliding_graph(["a", "b"]);
    graph.paths.get_mut(&PathId(0)).unwrap().nodes = vec![];

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Node 0 is not in the walk of its path 0"),
      "unexpected error: {err}"
    );
  }

  /// Node ids are unique across the graph, so a node claimed by two walks makes its own `path_id`
  /// meaningless, and would have the node reconstructed into two different genomes.
  #[rstest]
  fn test_validate_rejects_node_walked_by_two_paths() {
    let mut graph = colliding_graph(["a", "b"]);
    graph.paths.get_mut(&PathId(1)).unwrap().nodes = vec![NodeId(0), NodeId(1)];

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Node 0 is walked by both path 0 and path 1"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_validate_rejects_node_missing_from_the_alignments_of_its_block() {
    let mut graph = colliding_graph(["a", "b"]);
    graph
      .blocks
      .insert(BlockId(0), PangraphBlock::new(BlockId(0), "ACGTACGT", btreemap! {}));

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(err.contains("Node 0 not found in block 0"), "unexpected error: {err}");
  }

  /// Reconstruction rotates a genome by the offset of its first node. An offset past the end of the
  /// genome used to reach `rotate_right`, which panics rather than reporting.
  #[rstest]
  fn test_validate_rejects_node_position_beyond_the_length_of_its_path() {
    let mut graph = colliding_graph(["a", "b"]);
    let node = &graph.nodes[&NodeId(0)];
    graph.nodes.insert(
      NodeId(0),
      PangraphNode::new(
        Some(NodeId(0)),
        node.block_id(),
        node.path_id(),
        node.strand(),
        (100, 8),
      ),
    );

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Node 0 has position (100, 8), outside its path 0 of total length 8"),
      "unexpected error: {err}"
    );
  }

  /// `Edit::apply` indexes the consensus by edit position, so an out-of-range edit used to panic
  /// while reconstructing the block, even though `apply` returns a `Result`.
  #[rstest]
  fn test_validate_rejects_edit_beyond_the_consensus_of_its_block() {
    let mut graph = colliding_graph(["a", "b"]);
    let edit = Edit {
      subs: vec![Sub::new(99, 'T')],
      dels: vec![],
      inss: vec![],
    };
    graph.blocks.insert(
      BlockId(0),
      PangraphBlock::new(BlockId(0), "ACGTACGT", btreemap! { NodeId(0) => edit }),
    );

    let err = report_to_string(&graph.validate().unwrap_err());
    assert!(
      err.contains("Substitution position 99 is out of bounds for sequence of length 8"),
      "unexpected error: {err}"
    );
    assert!(
      err.contains("alignment of node 0 against block 0"),
      "unexpected error: {err}"
    );
  }
}
