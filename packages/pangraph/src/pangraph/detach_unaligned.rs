use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::representation::seq::Seq;
use crate::utils::id::id;
use eyre::Report;
use log::debug;
use std::collections::BTreeMap;

// takes as input a list of recently-merged blocks
// 1. looks for unaligned nodes
// 2. removes them from their parent blocks
// 3. creates new singleton blocks for each unaligned node
// 4. appends these new blocks to the list of blocks
// 5. returns a list of new node updates (to replace the old nodes in the pangraph)
pub fn detach_unaligned_nodes(
  blocks: &mut Vec<PangraphBlock>,
  nodes_dict: &mut BTreeMap<NodeId, PangraphNode>,
) -> Result<(), Report> {
  // Identify unaligned nodes and remove them from their blocks
  let unaligned_nodes = extract_unaligned_nodes(blocks)?;

  for (node_id, sequence) in unaligned_nodes {
    // Generate a unique block ID based on the node ID and sequence
    let new_block_id = BlockId(id((node_id, &sequence)));

    // Create a singleton block with the unaligned node
    let new_block = PangraphBlock::from_consensus(sequence, new_block_id, node_id);

    debug!(
      "Created new singleton block {} for unaligned node {}",
      new_block_id, node_id
    );
    // append the new block to the list of blocks
    blocks.push(new_block);

    // Create a new PangraphNode for the unaligned node
    let old_node = nodes_dict
      .get(&node_id)
      .ok_or_else(|| eyre::eyre!("Node {} not found in nodes dictionary", node_id))?;
    let new_node = PangraphNode::new(
      Some(node_id),
      old_node.block_id(),
      old_node.path_id(),
      old_node.strand(),
      old_node.position(),
    );

    // replace in the nodes dictionary
    nodes_dict.insert(node_id, new_node.clone());
    debug!("Replaced node {} in nodes dictionary with new node", node_id);
  }

  Ok(())
}

// Identify nodes that are unaligned, i.e. only consist of indels
// removes them from their blocks
// returns a list of node ids and their sequences
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
