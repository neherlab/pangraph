use crate::circularize::circularize_utils::{Edge, SimpleNode};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use eyre::Report;
use maplit::btreemap;
use std::borrow::Cow;
use std::collections::BTreeMap;

/// Merges the two blocks from the same edge into a single block,
/// modifying the graph in place.
pub fn merge_blocks(graph: &mut Pangraph, edge: Edge) -> Result<(), Report> {
  // 0) orient edge: deterministically assign a left and right block
  let edge = orient_merging_edge(graph, &edge);

  // 1) find node pairs: recovers pairs of nodes to merge,
  // and creates new nodes to be added to the graph
  let (node_map, new_nodes) = find_node_pairings(graph, &edge);

  // 2) merge alignment: concatenates the alignment of the two blocks
  // and returns a new block
  let new_block = merge_alignment(graph, &edge, &node_map, &new_nodes)?;

  // 3) update paths and nodes
  graph_merging_update(graph, new_block, &new_nodes, &edge);

  Ok(())
}

/// Returns the oriented edge for merging, with the first block being the anchor block.
/// This is:
/// - the block with the longest consensus sequence
/// - if the consensus lengths are equal, the block with the smallest id
fn orient_merging_edge(graph: &Pangraph, edge: &Edge) -> Edge {
  let b1 = &graph.blocks[&edge.n1.bid];
  let b2 = &graph.blocks[&edge.n2.bid];
  let l1 = b1.consensus_len();
  let l2 = b2.consensus_len();
  if l1 > l2 || (l1 == l2 && b1.id() < b2.id()) {
    *edge
  } else {
    edge.invert()
  }
}

/// Given a specific edge between blocks, returns pairings for all nodes in the blocks,
/// and a dictionary of new nodes to be added to the graph.
fn find_node_pairings(graph: &Pangraph, edge: &Edge) -> (BTreeMap<NodeId, NodeId>, BTreeMap<NodeId, PangraphNode>) {
  let mut node_pairings = btreemap! {};
  let mut new_nodes = btreemap! {};
  for (&path_id, path) in &graph.paths {
    let n = path.nodes.len();
    let i = if path.circular { n } else { n - 1 };
    for idx in 0..i {
      let nid1 = path.nodes[idx];
      let nid2 = path.nodes[(idx + 1) % n];
      let n1 = &graph.nodes[&nid1];
      let n2 = &graph.nodes[&nid2];

      let bid1 = n1.block_id();
      let bid2 = n2.block_id();
      let strand1 = n1.strand();
      let strand2 = n2.strand();

      let sn1 = SimpleNode::new(bid1, strand1);
      let sn2 = SimpleNode::new(bid2, strand2);
      let e = Edge { n1: sn1, n2: sn2 };

      if *edge == e {
        node_pairings.insert(nid1, nid2);
        node_pairings.insert(nid2, nid1);
        let (new_s, new_e) = (n1.position().0, n2.position().1);
        let new_strand = if edge.n1 == sn1 { n1.strand() } else { n2.strand() };
        debug_assert!(
          n1.position().1 % path.tot_len() == n2.position().0 % path.tot_len(),
          "nodes should be adjacent:\nn1: {:?},\nn2: {:?}",
          n1,
          n2
        );

        let new_node = PangraphNode::new(None, edge.n1.bid, path_id, new_strand, (new_s, new_e));
        new_nodes.insert(nid1, new_node.clone());
        new_nodes.insert(nid2, new_node);
      }
    }
  }
  (node_pairings, new_nodes)
}

/// Merges the alignment of the two blocks and returns a new alignment and consensus.
fn merge_alignment(
  graph: &Pangraph,
  edge: &Edge,
  node_map: &BTreeMap<NodeId, NodeId>,
  new_nodes: &BTreeMap<NodeId, PangraphNode>,
) -> Result<PangraphBlock, Report> {
  let new_node_ids: BTreeMap<NodeId, NodeId> = new_nodes.iter().map(|(&k, n)| (k, n.id())).collect();

  let b1 = &graph.blocks[&edge.n1.bid];
  let b2 = &graph.blocks[&edge.n2.bid];

  let (b_left, b_right) = if edge.n1.strand.is_forward() {
    (b1, b2)
  } else {
    (b2, b1)
  };

  let b_right = if edge.n1.strand != edge.n2.strand {
    Cow::Owned(b_right.reverse_complement()?)
  } else {
    Cow::Borrowed(b_right)
  };

  Ok(concatenate_alignments(
    b_left,
    &b_right,
    node_map,
    &new_node_ids,
    edge.n1.bid,
  ))
}

/// Concatenates two blocks, with new node ids.
fn concatenate_alignments(
  bl1: &PangraphBlock,
  bl2: &PangraphBlock,
  node_map: &BTreeMap<NodeId, NodeId>,
  new_nodes_id: &BTreeMap<NodeId, NodeId>,
  new_block_id: BlockId,
) -> PangraphBlock {
  debug_assert!(bl1.depth() == bl2.depth(), "blocks must have the same depth");

  let seq = [bl1.consensus(), bl2.consensus()].concat();

  let mut aln = btreemap! {};
  for (&nid1, e1) in bl1.alignments() {
    let nid2 = node_map[&nid1];
    let e2 = &bl2.alignment(nid2);
    let new_id = new_nodes_id[&nid1];
    aln.insert(new_id, e1.concat(&e2.shift(bl1.consensus_len() as isize)));
  }

  PangraphBlock::new(new_block_id, seq, aln)
}

/// Updates the graph in place after merging two blocks.
fn graph_merging_update(
  graph: &mut Pangraph,
  new_block: PangraphBlock,
  new_nodes: &BTreeMap<NodeId, PangraphNode>,
  edge: &Edge,
) {
  graph.blocks.remove(&edge.n1.bid);
  graph.blocks.remove(&edge.n2.bid);

  graph.blocks.insert(new_block.id(), new_block);

  let bid_left = edge.n1.bid;

  graph_merging_update_paths(graph, new_nodes, bid_left);

  graph_merging_update_nodes(graph, new_nodes, bid_left);
}

fn graph_merging_update_paths(graph: &mut Pangraph, new_nodes: &BTreeMap<NodeId, PangraphNode>, bid_left: BlockId) {
  for path in &mut graph.paths.values_mut() {
    path.nodes.retain_mut(|nid| match new_nodes.get(nid) {
      Some(new_node) if graph.nodes[nid].block_id() == bid_left => {
        *nid = new_node.id();
        true
      }
      Some(_) => false,
      None => true,
    });
  }
}

fn graph_merging_update_nodes(graph: &mut Pangraph, new_nodes: &BTreeMap<NodeId, PangraphNode>, bid_left: BlockId) {
  for (&nid, n) in new_nodes {
    if graph.nodes[&nid].block_id() == bid_left {
      graph.nodes.insert(n.id(), n.clone());
    }
    graph.nodes.remove(&nid);
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::circularize::circularize::remove_transitive_edges;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_path::{PangraphPath, PathId};
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  fn block_1() -> PangraphBlock {
    //          0         1         2         3
    //          01234567890123456789012345678901
    // cons:    ACTATATTACGGCGATCGATCGATTACTCGCT
    //   n1:    ...G............................  l = 32
    //   n2:    .......|.....xxx................  l = 31
    //   n3:    ................................| l = 35
    #[rustfmt::skip]
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![],                    vec![],                vec![Sub::new(3, 'G')]),
      NodeId(2) => Edit::new(vec![Ins::new(7, "AA")],   vec![Del::new(13, 3)], vec![]                ),
      NodeId(3) => Edit::new(vec![Ins::new(32, "CCC")], vec![],                vec![]                ),
    };
    PangraphBlock::new(BlockId(1), "ACTATATTACGGCGATCGATCGATTACTCGCT", aln)
  }

  fn block_2() -> PangraphBlock {
    //          0         1         2         3
    //          0123456789012345678901234567890
    // cons:    GATCTTAGGATCATCCCTATCATAGGAGTCG
    //   n4:    .........................xx....  l = 29
    //   n5:    ...T...........................  l = 31
    //   n6:   |xx.............................  l = 32
    #[rustfmt::skip]
    let aln = btreemap! {
      NodeId(4) => Edit::new(vec![],                   vec![Del::new(25, 2)],   vec![]                ),
      NodeId(5) => Edit::new(vec![],                   vec![],                  vec![Sub::new(3, 'T')]),
      NodeId(6) => Edit::new(vec![Ins::new(0, "TTT")], vec![Del::new(0, 2)],    vec![]                ),
    };
    PangraphBlock::new(BlockId(2), "GATCTTAGGATCATCCCTATCATAGGAGTCG", aln)
  }

  fn block_3() -> PangraphBlock {
    //          0         1         2
    //          012345678901234567890
    // cons:    CTATTACTAGGGGGACCACTA
    //   n7:    ...............xx....  l = 19
    //   n8:    ...C.................  l = 21
    #[rustfmt::skip]
    let aln = btreemap! {
      NodeId(7) => Edit::new(vec![],                   vec![Del::new(15, 2)],  vec![]                ),
      NodeId(8) => Edit::new(vec![],                   vec![],                 vec![Sub::new(3, 'C')]),
    };
    PangraphBlock::new(BlockId(3), "CTATTACTAGGGGGACCACTA", aln)
  }

  fn block_1_revcomp() -> PangraphBlock {
    //          0         1         2         3
    //          01234567890123456789012345678901
    // cons:    AGCGAGTAATCGATCGATCGCCGTAATATAGT
    //   n1:    ............................C...  l = 32
    //   n2:    ................xxx......|......  l = 31
    //   n3:    |...............................  l = 35
    PangraphBlock::new(
      BlockId(1),
      "AGCGAGTAATCGATCGATCGCCGTAATATAGT",
      btreemap! {
        NodeId(1) => Edit::new(vec![],                   vec![],                vec![Sub::new(28, 'C')]),
        NodeId(2) => Edit::new(vec![Ins::new(25, "TT")], vec![Del::new(16, 3)], vec![]                 ),
        NodeId(3) => Edit::new(vec![Ins::new(0, "GGG")], vec![],                vec![]                 ),
      },
    )
  }

  fn block_2_revcomp() -> PangraphBlock {
    //          0         1         2         3
    //          0123456789012345678901234567890
    // cons:    CGACTCCTATGATAGGGATGATCCTAAGATC
    //   n4:    ....xx.........................  l = 29
    //   n5:    ...........................A...  l = 31
    //   n6:    .............................xx| l = 32
    PangraphBlock::new(
      BlockId(2),
      "CGACTCCTATGATAGGGATGATCCTAAGATC",
      btreemap! {
        NodeId(4) => Edit::new(vec![],                    vec![Del::new(4, 2)],  vec![]                 ),
        NodeId(5) => Edit::new(vec![],                    vec![],                vec![Sub::new(27, 'A')]),
        NodeId(6) => Edit::new(vec![Ins::new(31, "AAA")], vec![Del::new(29, 2)], vec![]                 ),
      },
    )
  }

  fn graph_a() -> Pangraph {
    //      (0|32)      (32|61)     (61|0)
    // p1) (b1+|n1) -> (b2-|n4) -> (b3+|n7)  l=80
    //      (10|41)     (41|72)     (72|10)
    // p2) (b1+|n2) -> (b2-|n5) -> (b3+|n8)  l=83
    //      (5|40)      (40|5)
    // p3) (b2+|n6) -> (b1-|n3)              l=67

    #[rustfmt::skip]
    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), [1, 4, 7].into_iter().map(NodeId).collect_vec(), 80, true,  None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [2, 5, 8].into_iter().map(NodeId).collect_vec(), 83, true,  None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [6, 3].into_iter().map(NodeId).collect_vec(),    67, true, None),
    };

    let blocks = btreemap! {
      BlockId(1) => block_1(),
      BlockId(2) => block_2(),
      BlockId(3) => block_3(),
    };

    let nodes = btreemap! {
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(1), PathId(1), Forward, (0,  32)),
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), BlockId(1), PathId(2), Forward, (10, 41)),
      NodeId(3) => PangraphNode::new(Some(NodeId(3)), BlockId(1), PathId(3), Reverse, (40, 5 )),
      NodeId(4) => PangraphNode::new(Some(NodeId(4)), BlockId(2), PathId(1), Reverse, (32, 61)),
      NodeId(5) => PangraphNode::new(Some(NodeId(5)), BlockId(2), PathId(2), Reverse, (41, 72)),
      NodeId(6) => PangraphNode::new(Some(NodeId(6)), BlockId(2), PathId(3), Forward, (5,  40)),
      NodeId(7) => PangraphNode::new(Some(NodeId(7)), BlockId(3), PathId(1), Forward, (61, 0 )),
      NodeId(8) => PangraphNode::new(Some(NodeId(8)), BlockId(3), PathId(2), Forward, (72, 10)),
    };

    Pangraph { paths, blocks, nodes }
  }

  #[test]
  fn test_find_node_pairings_a() {
    let n1 = SimpleNode::new(BlockId(1), Forward);
    let n2 = SimpleNode::new(BlockId(2), Reverse);
    let edge = Edge { n1, n2 };

    let (pairings, new_nodes) = find_node_pairings(&graph_a(), &edge);

    assert_eq!(
      pairings,
      btreemap! {
        NodeId(1) => NodeId(4),
        NodeId(2) => NodeId(5),
        NodeId(3) => NodeId(6),
        NodeId(4) => NodeId(1),
        NodeId(5) => NodeId(2),
        NodeId(6) => NodeId(3),
      }
    );
  }

  #[test]
  fn test_reverse_complement_1() {
    assert_eq!(block_1_revcomp(), block_1().reverse_complement().unwrap());
  }

  #[test]
  fn test_reverse_complement_2() {
    assert_eq!(block_2_revcomp(), block_2().reverse_complement().unwrap());
  }

  fn expected_new_nodes() -> BTreeMap<NodeId, PangraphNode> {
    btreemap! {
      NodeId(1) => PangraphNode::new(None, BlockId(1), PathId(1), Forward, (0 , 61)),
      NodeId(2) => PangraphNode::new(None, BlockId(1), PathId(2), Forward, (10, 72)),
      NodeId(3) => PangraphNode::new(None, BlockId(1), PathId(3), Reverse, (5, 5)),
      NodeId(4) => PangraphNode::new(None, BlockId(1), PathId(1), Forward, (0 , 61)),
      NodeId(5) => PangraphNode::new(None, BlockId(1), PathId(2), Forward, (10, 72)),
      NodeId(6) => PangraphNode::new(None, BlockId(1), PathId(3), Reverse, (5, 5)),
    }
  }

  fn expected_new_node_ids() -> BTreeMap<NodeId, NodeId> {
    expected_new_nodes().iter().map(|(&k, v)| (k, v.id())).collect()
  }

  fn expected_concat() -> PangraphBlock {
    //          0         1         2         3         4         5         6
    //          012345678901234567890123456789012345678901234567890123456789012
    // cons:    ACTATATTACGGCGATCGATCGATTACTCGCTCGACTCCTATGATAGGGATGATCCTAAGATC
    //   n1:    ...G................................xx.........................  l = 32 + 29
    //   n2:    .......|.....xxx...........................................A...  l = 31 + 31
    //   n3:    ................................|............................xx| l = 35 + 32
    let new_ids = expected_new_node_ids();
    let aln = btreemap! {
      new_ids[&NodeId(1)] => Edit::new(vec![],                                         vec![Del::new(36, 2)], vec![Sub::new(3 , 'G')]),
      new_ids[&NodeId(2)] => Edit::new(vec![Ins::new(7, "AA")],                        vec![Del::new(13, 3)], vec![Sub::new(59, 'A')]),
      new_ids[&NodeId(3)] => Edit::new(vec![Ins::new(32, "CCC"), Ins::new(63, "AAA")], vec![Del::new(61, 2)], vec![]                 ),
    };

    PangraphBlock::new(
      BlockId(1),
      "ACTATATTACGGCGATCGATCGATTACTCGCTCGACTCCTATGATAGGGATGATCCTAAGATC",
      aln,
    )
  }

  #[test]
  fn test_concatenate_blocks() {
    let n1 = SimpleNode::new(BlockId(1), Forward);
    let n2 = SimpleNode::new(BlockId(2), Reverse);
    let edge = Edge::new(n1, n2);

    let (pairings, new_nodes) = find_node_pairings(&graph_a(), &edge);

    let new_nodes_ids = expected_new_node_ids();
    let block = concatenate_alignments(
      &block_1(),
      &block_2().reverse_complement().unwrap(),
      &pairings,
      &new_nodes_ids,
      BlockId(1),
    );
    assert_eq!(block, expected_concat());
  }

  fn expected_graph() -> Pangraph {
    //       (0|-----------|61)     (61|0)
    // p1) (b1+|-----------|n1) -> (b3+|n7)  l=80
    //      (10|-----------|72)     (72|10)
    // p2) (b1+|-----------|n2) -> (b3+|n8)  l=83
    //      (40|-----------|40)
    // p3) (b1-|-----------|n3)              l=67

    let new_ids = expected_new_node_ids();

    let blocks = btreemap! {
      BlockId(1) => expected_concat(),
      BlockId(3) => PangraphBlock::new(BlockId(3), "CTATTACTAGGGGGACCACTA", btreemap! {
        NodeId(7) => Edit::new(vec![], vec![Del::new(15, 2)], vec![]                ),
        NodeId(8) => Edit::new(vec![], vec![],                vec![Sub::new(3, 'C')]),
      })
    };

    let nodes = btreemap! {
      new_ids[&NodeId(1)] => PangraphNode::new(None, BlockId(1), PathId(1), Forward, (0 , 61)),
      new_ids[&NodeId(2)] => PangraphNode::new(None, BlockId(1), PathId(2), Forward, (10, 72)),
      new_ids[&NodeId(3)] => PangraphNode::new(None, BlockId(1), PathId(3), Reverse, (5, 5)),
      NodeId(7) => PangraphNode::new(Some(NodeId(7)), BlockId(3), PathId(1), Forward, (61, 0 )),
      NodeId(8) => PangraphNode::new(Some(NodeId(8)), BlockId(3), PathId(2), Forward, (72, 10)),
    };

    #[rustfmt::skip]
    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), [new_ids[&NodeId(1)], NodeId(7)], 80, true, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [new_ids[&NodeId(2)], NodeId(8)], 83, true, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [new_ids[&NodeId(3)]],            67, true, None)
    };

    Pangraph { paths, blocks, nodes }
  }

  #[test]
  fn test_update_paths() {
    let n1 = SimpleNode::new(BlockId(1), Forward);
    let n2 = SimpleNode::new(BlockId(2), Reverse);
    let edge = Edge::new(n1, n2);
    let mut g = graph_a();
    let (pairings, new_nodes) = find_node_pairings(&g, &edge);
    graph_merging_update_paths(&mut g, &new_nodes, BlockId(1));
    assert_eq!(g.paths, expected_graph().paths);
  }

  #[test]
  fn test_update_nodes() {
    let n1 = SimpleNode::new(BlockId(1), Forward);
    let n2 = SimpleNode::new(BlockId(2), Reverse);
    let edge = Edge::new(n1, n2);
    let mut g = graph_a();
    let (pairings, new_nodes) = find_node_pairings(&g, &edge);
    graph_merging_update_nodes(&mut g, &new_nodes, BlockId(1));
    assert_eq!(g.nodes, expected_graph().nodes);
  }

  #[test]
  fn test_merge_blocks() {
    let n1 = SimpleNode::new(BlockId(1), Forward);
    let n2 = SimpleNode::new(BlockId(2), Reverse);
    let edge = Edge::new(n1, n2);
    let mut g = graph_a();
    merge_blocks(&mut g, edge).unwrap();
    assert_eq!(g.nodes, expected_graph().nodes);
    assert_eq!(g.blocks, expected_graph().blocks);
    assert_eq!(g.paths, expected_graph().paths);
  }

  #[test]
  fn test_remove_transitive_edges() {
    let mut g = graph_a();
    remove_transitive_edges(&mut g).unwrap();
    assert_eq!(g.nodes, expected_graph().nodes);
    assert_eq!(g.blocks, expected_graph().blocks);
    assert_eq!(g.paths, expected_graph().paths);
  }
}
