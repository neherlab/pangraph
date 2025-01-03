use log::info;

use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_node::NodeId;

/// If present, removes empty nodes from a graph inplace. Only specified
/// blocks are checked for empty nodes.
///
/// Remove empty nodes:
/// - for blocks that originate from a new merger or split
/// - goes through each block and removes any node with an empty sequence
/// - also remove the node from the graph node list and paths
pub fn remove_emtpy_nodes(graph: &mut Pangraph, block_ids: &[BlockId]) {
  let node_ids = find_empty_nodes(graph, block_ids);
  remove_nodes_from_graph(graph, &node_ids);
}

/// Finds nodes with empty sequences (full deletion) in the graph.
/// It only checks specific blocks (the ones that were updated by a merger/split).
fn find_empty_nodes(graph: &Pangraph, block_ids: &[BlockId]) -> Vec<NodeId> {
  let mut node_ids_to_delete = Vec::new();
  for &block_id in block_ids {
    let block = &graph.blocks[&block_id];
    let cons_len = block.consensus_len();

    debug_assert!(
      cons_len > 0,
      "Block {block_id} has a consensus length of 0 and should have been removed",
    );

    for (&node_id, edits) in block.alignments() {
      // check that edits are non-overlapping and well-defined
      #[cfg(any(test, debug_assertions))]
      edits.sanity_check(cons_len).unwrap();

      //  if the node has insertions or substitutions, or the total deletion size is less than the consensus length
      // then it is not empty
      if !edits.inss.is_empty() || !edits.subs.is_empty() || edits.dels.is_empty() {
        // check that the node is not empty
        // first exclude edge-case: circular path with a single node. In this case start == end but the node is not empty
        let path_id = graph.nodes[&node_id].path_id();
        let path = &graph.paths[&path_id];
        let is_circular = path.circular();
        let single_node = path.nodes.len() == 1;
        if is_circular && single_node {
          debug_assert!(
            path.tot_len > 0,
            "Circular path {path_id} with a single node should not have a length of 0"
          );
          continue;
        }

        // if this is not the case, then it is sufficient to check that the node start != end
        debug_assert!(
          !graph.nodes[&node_id].start_is_end(),
          "Node {node_id} with edits {edits:?} and consensus length {cons_len} is empty and should have been removed",
        );

        continue;
      }

      // if the node has no insertions or substitutions and the total deletion size is equal to the consensus length
      // then it is empty and should be removed
      if edits.dels.iter().map(|d| d.len).sum::<usize>() == cons_len {
        info!(
          "empty node {} with edits {:?} and consensus length {} is empty: to be removed",
          node_id, edits, cons_len
        );
        debug_assert!(graph.nodes[&node_id].start_is_end());
        node_ids_to_delete.push(node_id);
      }
    }
  }
  node_ids_to_delete
}

/// Removes each node from the graph inplace. The node gets removed from:
/// - the node dictionary
/// - the alignment object in block
/// - the path
fn remove_nodes_from_graph(graph: &mut Pangraph, node_ids: &[NodeId]) {
  for &node_id in node_ids {
    let node = &graph.nodes[&node_id];
    let path_id = node.path_id();
    let block_id = node.block_id();

    // remove from node dictionary
    graph.nodes.remove(&node_id);

    // === old implementation. Assumes only one node to remove per path ===
    // === but there currently is an edge case where multiple empty nodes with the same id could be present ===
    // // remove from path
    // let path_nodes = &mut graph.paths.get_mut(&path_id).unwrap().nodes;
    // let node_idx = path_nodes.iter().position(|&n| n == node_id).unwrap();
    // path_nodes.remove(node_idx);

    //  remove all nodes with the same id
    let path_nodes = &mut graph.paths.get_mut(&path_id).unwrap().nodes;
    // while there is a node with this id, remove it:
    while let Some(node_idx) = path_nodes.iter().position(|&n| n == node_id) {
      path_nodes.remove(node_idx);
    }

    // remove from block alignment
    let _removed = graph.blocks.get_mut(&block_id).unwrap().alignment_remove(node_id);
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Sub};
  use crate::pangraph::pangraph_block::PangraphBlock;
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::pangraph_path::{PangraphPath, PathId};
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::pretty_assert_eq;
  use maplit::btreemap;

  fn create_input_graph() -> Pangraph {
    let nodes = btreemap! {
        NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(1), PathId(0), Forward, (0, 10)),
        NodeId(2) => PangraphNode::new(Some(NodeId(2)), BlockId(1), PathId(1), Forward, (0, 10)),
        NodeId(3) => PangraphNode::new(Some(NodeId(3)), BlockId(1), PathId(2), Reverse, (0, 0)),
        NodeId(4) => PangraphNode::new(Some(NodeId(4)), BlockId(2), PathId(0), Forward, (10, 20)),
        NodeId(5) => PangraphNode::new(Some(NodeId(5)), BlockId(2), PathId(2), Forward, (0, 10)),
    };

    let paths = btreemap! {
        PathId(0) => PangraphPath::new(Some(PathId(1)), vec![NodeId(1), NodeId(4)], 20, false, None),
        PathId(1) => PangraphPath::new(Some(PathId(2)), vec![NodeId(2)],            10, false, None),
        PathId(2) => PangraphPath::new(Some(PathId(3)), vec![NodeId(3), NodeId(5)], 10, false, None),
    };

    let blocks = btreemap! {
        BlockId(1) => PangraphBlock::new(BlockId(1), "AAAAAAAAAA", btreemap!{
          NodeId(1) => Edit::new(vec![], vec![Del::new(1, 3)],  vec![]),
          NodeId(2) => Edit::new(vec![], vec![],                vec![Sub::new(5, 'G')]),
          NodeId(3) => Edit::new(vec![], vec![Del::new(0, 10)], vec![]),
        }),
        BlockId(2) => PangraphBlock::new(BlockId(2), "CCCCCCCCCC", btreemap!{
          NodeId(4) => Edit::empty(),
          NodeId(5) => Edit::empty(),
        }),
    };

    Pangraph { paths, blocks, nodes }
  }

  #[test]
  fn test_find_empty_nodes() {
    let graph = create_input_graph();
    let empty_node_ids = find_empty_nodes(&graph, &[BlockId(1), BlockId(2)]);
    pretty_assert_eq!(empty_node_ids, vec![NodeId(3)]);
  }

  #[test]
  fn test_remove_empty_nodes() {
    let mut graph = create_input_graph();
    remove_emtpy_nodes(&mut graph, &[BlockId(1), BlockId(2)]);
    pretty_assert_eq!(graph.nodes.keys(), &[NodeId(1), NodeId(2), NodeId(4), NodeId(5)]);
    pretty_assert_eq!(graph.paths[&PathId(2)].nodes, vec![NodeId(5)]);
    assert!(!graph.blocks[&BlockId(1)].alignments().contains_key(&NodeId(3)));
  }
}
