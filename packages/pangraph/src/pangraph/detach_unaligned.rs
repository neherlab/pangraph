use crate::io::seq::reverse_complement;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::pangraph::strand::Strand::Forward;
use crate::representation::seq::Seq;
use crate::utils::id::id;
use eyre::Report;
use log::debug;
use std::collections::BTreeMap;

/// removes unaligned nodes from the graph. These are nodes whose alignments consist only of indels.
/// takes as input a list of recently-merged blocks and the nodes dictionary.
/// 1. looks for unaligned nodes
/// 2. removes them from their parent blocks
/// 3. creates new singleton blocks for each unaligned node
/// 4. appends these new blocks to the list of blocks
/// 5. substitute the new nodes in the nodes dictionary, keeping the same node_id so that paths are not affected
pub fn detach_unaligned_nodes(
  blocks: &mut Vec<PangraphBlock>,
  nodes_dict: &mut BTreeMap<NodeId, PangraphNode>,
) -> Result<(), Report> {
  // Identify unaligned nodes and remove them from their blocks
  let unaligned_nodes = extract_unaligned_nodes(blocks)?;

  for (node_id, sequence) in unaligned_nodes {
    // removes the old node from the nodes dictionary
    let old_node = nodes_dict
      .remove(&node_id)
      .ok_or_else(|| eyre::eyre!("Node {} not found in nodes dictionary", node_id))?;
    // create new node and block
    let (new_node, new_block) = create_new_node_and_block(node_id, sequence, &old_node)?;
    // add the new block to the blocks list
    blocks.push(new_block);
    // replace the new node in the nodes dictionary
    nodes_dict.insert(node_id, new_node.clone());
    debug!("Replaced node {} in nodes dictionary with new node", node_id);
  }

  Ok(())
}

/// Identify nodes that are unaligned, i.e. only consist of indels
/// removes them from their blocks
/// returns a list of node ids and their sequences
fn extract_unaligned_nodes(blocks: &mut [PangraphBlock]) -> Result<Vec<(NodeId, Seq)>, Report> {
  let mut unaligned_nodes = Vec::new();

  for block in blocks.iter_mut() {
    let cons_len = block.consensus_len();

    // Find nodes that are unaligned (only consist of indels)
    let mut removed = Vec::new();
    for (&node_id, edit) in block.alignments() {
      if edit.aligned_count(cons_len) == 0 {
        debug!("Found unaligned node {} in block {}", node_id, block.id());
        removed.push(node_id);
      }
    }
    // Remove unaligned nodes from the block and collect them
    for node_id in removed {
      let edit = block.alignment_remove(node_id);
      let seq = edit.apply(block.consensus())?;
      unaligned_nodes.push((node_id, seq));
    }
  }
  Ok(unaligned_nodes)
}

/// Create a new PangraphNode and PangraphBlock for an unaligned node
/// Uses the same node_id, and optionally reverse-complements the sequence
/// if the original node was on the reverse strand.
fn create_new_node_and_block(
  node_id: NodeId,
  seq: Seq,
  old_node: &PangraphNode,
) -> Result<(PangraphNode, PangraphBlock), Report> {
  // if the node is on the reverse strand, we need to reverse the sequence
  let seq = if old_node.strand() == Forward {
    seq
  } else {
    reverse_complement(&seq)?
  };

  // Generate a unique block ID based on the node ID and sequence
  let new_block_id = BlockId(id((node_id, &seq)));

  // Create a singleton block with the unaligned node
  let new_block = PangraphBlock::from_consensus(seq, new_block_id, node_id);

  // Create a new PangraphNode for the unaligned node
  let new_node = PangraphNode::new(
    Some(node_id),
    new_block_id,
    old_node.path_id(),  // same path ID as the old node
    Forward,             // assuming the new node is always on the forward strand
    old_node.position(), // same position as the old node
  );

  Ok((new_node, new_block))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
  use crate::pangraph::pangraph_path::PathId;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::representation::seq::Seq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_create_new_node_and_block_forward() -> Result<(), Report> {
    let node_id = NodeId(1);
    let seq = Seq::from_str("ATGTTGATAG");
    let old_block_id = BlockId(0);
    let old_path_id = PathId(0);
    let old_node = PangraphNode::new(Some(node_id), old_block_id, old_path_id, Forward, (10, 20));

    let (new_node, new_block) = create_new_node_and_block(node_id, seq.clone(), &old_node)?;

    let expected_new_node = PangraphNode::new(Some(node_id), new_block.id(), old_path_id, Forward, (10, 20));
    let expected_new_block = PangraphBlock::from_consensus(seq, new_block.id(), node_id);

    assert_eq!(new_node, expected_new_node);
    assert_eq!(new_block, expected_new_block);

    Ok(())
  }

  #[rstest]
  fn test_create_new_node_and_block_reverse() -> Result<(), Report> {
    let node_id = NodeId(2);
    let seq = Seq::from_str("ATGTTGATAG");
    let old_block_id = BlockId(0);
    let old_path_id = PathId(1);
    let old_node = PangraphNode::new(Some(node_id), old_block_id, old_path_id, Reverse, (5, 15));

    let (new_node, new_block) = create_new_node_and_block(node_id, seq.clone(), &old_node)?;

    let expected_new_node = PangraphNode::new(Some(node_id), new_block.id(), old_path_id, Forward, (5, 15));
    let expected_new_block = PangraphBlock::from_consensus(reverse_complement(&seq)?, new_block.id(), node_id);

    assert_eq!(new_node, expected_new_node);
    assert_eq!(new_block, expected_new_block);

    Ok(())
  }

  #[rstest]
  fn test_extract_unaligned_nodes_simple() -> Result<(), Report> {
    let block_id = BlockId(0);
    let block_cons = Seq::from_str("AAAAAAAAAAAAAAAA");
    // two nodes: one aligned, one unaligned
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![],                    vec![],                vec![Sub::new(3, 'G')]),
      NodeId(2) => Edit::new(vec![Ins::new(16, "CCCCCCCC")],   vec![Del::new(0, 16)], vec![]),
    };
    let block = PangraphBlock::new(block_id, block_cons.clone(), aln);

    let mut blocks = vec![block];
    let unaligned = extract_unaligned_nodes(blocks.as_mut_slice())?;

    // check that node_2 was extracted
    assert_eq!(unaligned, vec![(NodeId(2), Seq::from_str("CCCCCCCC"))]);
    // check that the node was removed from the block alignment
    assert_eq!(blocks.len(), 1);
    let expected_bl = PangraphBlock::new(
      block_id,
      block_cons,
      btreemap! {NodeId(1) => Edit::new(vec![], vec![], vec![Sub::new(3, 'G')])},
    );
    assert_eq!(blocks[0], expected_bl);

    Ok(())
  }

  #[rstest]
  fn test_detach_unaligned_nodes() -> Result<(), Report> {
    // setup one block with one aligned and one unaligned node
    let cons = Seq::from_str("AAAAAAAAAAAAAAAA");
    let aligned_edit = Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]);
    let unaligned_edit = Edit::new(vec![Ins::new(0, "CCCCCCCC")], vec![Del::new(0, 16)], vec![]);
    let block = PangraphBlock::new(
      BlockId(0),
      cons.clone(),
      btreemap! {
        NodeId(1) => aligned_edit.clone(),
        NodeId(2) => unaligned_edit,
      },
    );
    let mut blocks = vec![block];
    let node1 = PangraphNode::new(Some(NodeId(1)), BlockId(0), PathId(0), Forward, (0, 16));
    let node2 = PangraphNode::new(Some(NodeId(2)), BlockId(0), PathId(1), Reverse, (0, 8));
    let mut nodes = btreemap! {
      NodeId(1) => node1,
      NodeId(2) => node2,
    };
    // perform detach
    detach_unaligned_nodes(&mut blocks, &mut nodes)?;
    // original block now only has the aligned node
    assert_eq!(blocks.len(), 2);

    let expected_block1 = PangraphBlock::new(BlockId(0), cons, btreemap! {NodeId(1) => aligned_edit});
    let new_block_seq = Seq::from_str("GGGGGGGG");
    let new_block_id = BlockId(id((NodeId(2), &new_block_seq)));
    let expected_block2 = PangraphBlock::from_consensus(Seq::from_str("GGGGGGGG"), new_block_id, NodeId(2));

    assert_eq!(blocks[0], expected_block1);
    assert_eq!(blocks[1], expected_block2);

    let expected_node_dict = btreemap! {
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(0), PathId(0), Forward, (0, 16)),
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), new_block_id, PathId(1), Forward, (0, 8)),
    };
    assert_eq!(nodes, expected_node_dict);
    Ok(())
  }
}
