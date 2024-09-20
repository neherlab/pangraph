use crate::circularize::circularize_utils::{Edge, SimpleNode, SimplePath};
use crate::circularize::merge_blocks::merge_blocks;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_path::PathId;
use eyre::Report;
use std::collections::{BTreeMap, HashMap};

/// Removes transitive edges from the graph inplace.
pub fn remove_transitive_edges(graph: &mut Pangraph) -> Result<(), Report> {
  while let Some(edge) = find_transitive_edges(graph).first() {
    merge_blocks(graph, *edge)?;
  }
  Ok(())
}

/// Find transitive edges between two different blocks in the graph (no self-loops).
///
/// TODO: explain what transitive edges are
fn find_transitive_edges(graph: &Pangraph) -> Vec<Edge> {
  let block_depths = calculate_block_depths(graph);
  let edge_counts = count_edges(graph);
  let mut transitive_edges = Vec::new();

  #[allow(clippy::iter_over_hash_type)]
  for (edge, edge_count) in edge_counts {
    let bid1 = edge.n1.bid;
    let bid2 = edge.n2.bid;
    if block_depths[&bid1] == edge_count && block_depths[&bid2] == edge_count && bid1 != bid2 {
      transitive_edges.push(edge);
    }
  }
  transitive_edges
}

/// Map blocks ids to the block depths
fn calculate_block_depths(graph: &Pangraph) -> BTreeMap<BlockId, usize> {
  graph.blocks.iter().map(|(&bid, b)| (bid, b.depth())).collect()
}

/// Map edges to number of their occurences in the graph
fn count_edges(graph: &Pangraph) -> HashMap<Edge, usize> {
  let paths = graph_to_simple_paths(graph);
  let mut edges_ct = HashMap::new();
  for simple_path in paths.values() {
    for edge in simple_path.to_edges() {
      *edges_ct.entry(edge).or_insert(0) += 1;
    }
  }
  edges_ct
}

/// Map graph path ids to simple paths
fn graph_to_simple_paths(graph: &Pangraph) -> BTreeMap<PathId, SimplePath> {
  graph
    .paths
    .iter()
    .map(|(&path_id, path)| {
      let nodes = path
        .nodes
        .iter()
        .map(|node_id| {
          let node = &graph.nodes[node_id];
          SimpleNode::from_full_node(node)
        })
        .collect();

      (path_id, SimplePath::new(nodes, path.circular))
    })
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::PangraphBlock;
  use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
  use crate::pangraph::pangraph_path::PangraphPath;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  fn mock_aln(node_ids: impl IntoIterator<Item = usize>) -> BTreeMap<NodeId, Edit> {
    node_ids.into_iter().map(|nid| (NodeId(nid), Edit::empty())).collect()
  }

  fn input_graph() -> Pangraph {
    // create the following graph with circular paths:
    // a) 1+(10) -> 2+(20) -> 3+(30) -> 4+(40)
    // b) 1+(11) -> 2-(21) -> 2+(22) -> 3+(31) -> 4+(41)
    // c) 1+(12) -> 2+(23) -> 3-(32) -> 4+ (42)
    // d) 1+(13) -> 3-(33) -> 2+(24) -> 3-(34) -> 4+(43)
    // f) 4- (44) -> 3-(35) -> 2-(25) -> 1-(14)

    #[rustfmt::skip]
    let paths = btreemap!{
      PathId(0) => PangraphPath::new(Some(PathId(0)), [10, 20, 30, 40].into_iter().map(NodeId).collect_vec(),     0, true, None),
      PathId(1) => PangraphPath::new(Some(PathId(1)), [11, 21, 22, 31, 41].into_iter().map(NodeId).collect_vec(), 0, true, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [12, 23, 32, 42].into_iter().map(NodeId).collect_vec(),     0, true, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [13, 33, 24, 34, 43].into_iter().map(NodeId).collect_vec(), 0, true, None),
      PathId(4) => PangraphPath::new(Some(PathId(4)), [44, 35, 25, 14].into_iter().map(NodeId).collect_vec(),     0, true, None),
    };

    #[rustfmt::skip]
    let nodes = btreemap! {
      NodeId(10) => PangraphNode::new(Some(NodeId(10)), BlockId(1), PathId(0), Forward, (0, 0)),
      NodeId(20) => PangraphNode::new(Some(NodeId(20)), BlockId(2), PathId(0), Forward, (0, 0)),
      NodeId(30) => PangraphNode::new(Some(NodeId(30)), BlockId(3), PathId(0), Forward, (0, 0)),
      NodeId(40) => PangraphNode::new(Some(NodeId(40)), BlockId(4), PathId(0), Forward, (0, 0)),
      NodeId(11) => PangraphNode::new(Some(NodeId(11)), BlockId(1), PathId(1), Forward, (0, 0)),
      NodeId(21) => PangraphNode::new(Some(NodeId(21)), BlockId(2), PathId(1), Reverse, (0, 0)),
      NodeId(22) => PangraphNode::new(Some(NodeId(22)), BlockId(2), PathId(1), Forward, (0, 0)),
      NodeId(31) => PangraphNode::new(Some(NodeId(31)), BlockId(3), PathId(1), Forward, (0, 0)),
      NodeId(41) => PangraphNode::new(Some(NodeId(41)), BlockId(4), PathId(1), Forward, (0, 0)),
      NodeId(12) => PangraphNode::new(Some(NodeId(12)), BlockId(1), PathId(2), Forward, (0, 0)),
      NodeId(23) => PangraphNode::new(Some(NodeId(23)), BlockId(2), PathId(2), Forward, (0, 0)),
      NodeId(32) => PangraphNode::new(Some(NodeId(32)), BlockId(3), PathId(2), Reverse, (0, 0)),
      NodeId(42) => PangraphNode::new(Some(NodeId(42)), BlockId(4), PathId(2), Forward, (0, 0)),
      NodeId(13) => PangraphNode::new(Some(NodeId(13)), BlockId(1), PathId(3), Forward, (0, 0)),
      NodeId(33) => PangraphNode::new(Some(NodeId(33)), BlockId(3), PathId(3), Reverse, (0, 0)),
      NodeId(24) => PangraphNode::new(Some(NodeId(24)), BlockId(2), PathId(3), Forward, (0, 0)),
      NodeId(34) => PangraphNode::new(Some(NodeId(34)), BlockId(3), PathId(3), Reverse, (0, 0)),
      NodeId(43) => PangraphNode::new(Some(NodeId(43)), BlockId(4), PathId(3), Forward, (0, 0)),
      NodeId(44) => PangraphNode::new(Some(NodeId(44)), BlockId(4), PathId(4), Reverse, (0, 0)),
      NodeId(35) => PangraphNode::new(Some(NodeId(35)), BlockId(3), PathId(4), Reverse, (0, 0)),
      NodeId(25) => PangraphNode::new(Some(NodeId(25)), BlockId(2), PathId(4), Reverse, (0, 0)),
      NodeId(14) => PangraphNode::new(Some(NodeId(14)), BlockId(1), PathId(4), Reverse, (0, 0)),
    };

    #[rustfmt::skip]
    let blocks = btreemap!{
      BlockId(1) => PangraphBlock::new(BlockId(1), "", mock_aln([10, 11, 12, 13, 14])),
      BlockId(2) => PangraphBlock::new(BlockId(2), "", mock_aln([20, 21, 22, 23, 24, 25])),
      BlockId(3) => PangraphBlock::new(BlockId(3), "", mock_aln([30, 31, 32, 33, 34, 35])),
      BlockId(4) => PangraphBlock::new(BlockId(4), "", mock_aln([40, 41, 42, 43, 44])),
    };

    Pangraph { paths, blocks, nodes }
  }

  #[test]
  fn test_block_depths() {
    let actual = calculate_block_depths(&input_graph());
    let expected = btreemap! {BlockId(1) => 5, BlockId(2) => 6, BlockId(3) => 6, BlockId(4) => 5};
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_count_edges() {
    let ec = count_edges(&input_graph());

    let n1 = SimpleNode::new(BlockId(1), Forward);
    let n2 = SimpleNode::new(BlockId(2), Forward);
    let n3 = SimpleNode::new(BlockId(3), Forward);
    let n4 = SimpleNode::new(BlockId(4), Forward);

    assert_eq!(ec[&Edge::new(n1, n2)], 3);
    assert_eq!(ec[&Edge::new(n1, n2.invert())], 1);
    assert_eq!(ec[&Edge::new(n2, n3)], 3);
    assert_eq!(ec[&Edge::new(n2, n3.invert())], 2);
    assert_eq!(ec[&Edge::new(n2.invert(), n2)], 1);

    assert!(!ec.contains_key(&Edge::new(n2, n2.invert())));

    assert_eq!(ec[&Edge::new(n3, n4)], 3);
    assert_eq!(ec[&Edge::new(n3.invert(), n4)], 2);
    assert_eq!(ec[&Edge::new(n4, n1)], 5);
  }

  #[test]
  fn test_transitive_edges_a() {
    let n1 = SimpleNode::new(BlockId(1), Forward);
    let n4 = SimpleNode::new(BlockId(4), Forward);
    let expected = vec![Edge::new(n4, n1)];

    let actual = find_transitive_edges(&input_graph());
    assert_eq!(expected, actual);
  }

  fn block_1() -> PangraphBlock {
    //          0         1         2         3
    //          01234567890123456789012345678901
    // cons:    ACTATATTACGGCGATCGATCGATTACTCGCT
    //   n1:    ...G............................  l = 32
    //   n2:    .......|.....xxx................  l = 31
    //   n3:    ................................| l = 35
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![],                    vec![],                vec![Sub::new(3, 'G')]),
      NodeId(2) => Edit::new(vec![Ins::new(7, "AA")],   vec![Del::new(13, 3)], vec![]),
      NodeId(3) => Edit::new(vec![Ins::new(31, "CCC")], vec![],                vec![]),
    };
    PangraphBlock::new(BlockId(1), "ACTATATTACGGCGATCGATCGATTACTCGCT", aln)
  }

  fn graph_b() -> Pangraph {
    // single block: no mergings
    //      (0|32)
    // p1) (b1+|n1)  l=32
    //      (0|31)
    // p2) (b1+|n2)  l=31
    //      (0|35)
    // p3) (b1-|n3)  l=35
    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), vec![NodeId(1)], 32, true, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), vec![NodeId(2)], 31, true, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), vec![NodeId(3)], 35, true, None)
    };
    let blocks = btreemap! {
      BlockId(1) => block_1()
    };
    let nodes = btreemap! {
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(1), PathId(1), Forward, (0, 32)),
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), BlockId(1), PathId(2), Forward, (0, 31)),
      NodeId(3) => PangraphNode::new(Some(NodeId(3)), BlockId(1), PathId(3), Reverse, (0, 35))
    };
    Pangraph { paths, blocks, nodes }
  }

  #[test]
  fn test_transitive_edges_b() {
    let te = find_transitive_edges(&graph_b());
    assert_eq!(te, vec![]);
  }
}
