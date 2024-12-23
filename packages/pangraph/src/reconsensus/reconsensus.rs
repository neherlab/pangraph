use crate::align::map_variations::map_variations;
use crate::make_error;
use crate::pangraph::edits::{Ins, Sub};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::NodeId;
use crate::reconsensus::remove_nodes::remove_emtpy_nodes;
use crate::utils::collections::insert_at_inplace;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::BTreeMap;

/// Applies the reconsensus operation to each updated block in the graph:
/// - updates the block consensus following a merge
/// - removes potentially empty nodes
///
/// Reconsensus function:
///   - for blocks that:
///     - originate from a new merger
///     - have depth > 2
///   - goes through each block and re-defines the consensus
///     - for mutations, changes any mutated position. Also update the alignment (this is straightforward)
///     - for any in/del present in > M/2 sites, appends it to the consensus
///   - if the consensus has updated indels, then re-aligns all the sequences to the new consensus
pub fn reconsensus_graph(graph: &mut Pangraph, ids_updated_blocks: Vec<BlockId>) -> Result<(), Report> {
  // remove selected nodes from graph
  remove_emtpy_nodes(graph, &ids_updated_blocks);

  for block_id in ids_updated_blocks {
    let block = graph.blocks.get_mut(&block_id).unwrap();
    reconsensus(block)?;
  }

  Ok(())
}

/// Performs the reconsensus operation inplace on a block.
/// - if a position is mutated in > N/2 sites, it adds the mutation to the consensus and updates the alignment.
/// - if an in/del is present in > N/2 sites, it adds it to the consensus and re-aligns the sequences to the updated consensus.
fn reconsensus(block: &mut PangraphBlock) -> Result<(), Report> {
  reconsensus_mutations(block)?;
  let ins = majority_insertions(block);
  let dels = majority_deletions(block);
  if !ins.is_empty() || !dels.is_empty() {
    let consensus = block.consensus();
    let consensus = apply_indels(consensus, &dels, &ins);

    // debug assert: consensus is not empty
    debug_assert!(!consensus.is_empty(), "Consensus is empty after indels");

    update_block_consensus(block, &consensus)?;
  }
  Ok(())
}

/// Re-computes the consensus for a block if a position is mutated in > N/2 sites.
fn reconsensus_mutations(block: &mut PangraphBlock) -> Result<(), Report> {
  let n = block.depth();
  let mut muts = btreemap! {};

  // count mutations
  for (&nid, edit) in block.alignments() {
    for s in &edit.subs {
      *muts
        .entry(s.pos)
        .or_insert_with(BTreeMap::new)
        .entry(s.alt)
        .or_insert(0) += 1;
    }
  }

  // change positions that are different in more than N/2 sites.
  let mut changes = Vec::new();
  for (pos, ct) in muts {
    let (alt, count) = ct.iter().max_by_key(|&(_, v)| v).unwrap();
    if *count > n / 2 {
      changes.push((pos, *alt));
    }
  }

  // apply change
  for (pos, alt) in changes {
    let original = block.consensus().chars().nth(pos).unwrap();

    // update consensus
    block.set_consensus_char(pos, alt);

    // change mutations
    for edit in &mut block.alignments_mut().values_mut() {
      let subs_at_pos: Vec<_> = edit.subs.iter().filter(|s| s.pos == pos).cloned().collect();
      match subs_at_pos.len() {
        0 => {
          // if the position is not in a deletion, append the mutation
          if !edit.dels.iter().any(|d| d.contains(pos)) {
            edit.subs.push(Sub::new(pos, original));
            edit.subs.sort_by_key(|s| s.pos);
          }
        }
        1 => {
          let s = &subs_at_pos[0];
          if s.alt == alt {
            edit.subs.retain(|sub| !(sub.pos == pos && sub.alt == alt));
          }
        }
        _ => {
          return make_error!(
            "At block {}: at position {pos}: sequence states disagree: {:}",
            block.id(),
            subs_at_pos
              .iter()
              .map(|sub| sub.alt.to_string())
              .collect_vec()
              .join(", ")
          )
        }
      }
    }
  }

  Ok(())
}

/// Returns a list of positions to be removed from the consensus, because they are deleted in > N/2 sites.
fn majority_deletions(block: &PangraphBlock) -> Vec<usize> {
  let n = block.depth();
  let mut n_dels = btreemap! {};
  // for each deleted position, increment the counter
  for edit in block.alignments().values() {
    for d in &edit.dels {
      for pos in d.range() {
        *n_dels.entry(pos).or_insert(0) += 1;
      }
    }
  }
  // return the positions that are deleted in more than N/2 sites
  n_dels
    .into_iter()
    .filter(|&(_, count)| count > n / 2)
    .map(|(pos, _)| pos)
    .collect()
}

/// Returns a list of insertions to be added to the consensus, because they are inserted in > N/2 sites.
fn majority_insertions(block: &PangraphBlock) -> Vec<Ins> {
  let n = block.depth();
  let mut n_ins = btreemap! {};
  for edit in block.alignments().values() {
    for i in &edit.inss {
      *n_ins.entry((i.pos, i.seq.clone())).or_insert(0) += 1;
    }
  }

  // return the positions that are inserted in more than N/2 sites
  n_ins
    .into_iter()
    .filter(|&(_, count)| count > n / 2)
    .map(|((pos, ins), _)| Ins::new(pos, ins))
    .collect()
}

/// Updates the consensus sequence with the deletions and insertions.
fn apply_indels(cons: &str, dels: &[usize], inss: &[Ins]) -> String {
  let mut cons: Vec<char> = cons.chars().collect();

  for &pos in dels {
    cons[pos] = '\0'; // Using '\0' to temporarily denote deleted positions
  }

  // Reverse to maintain correct insertion indexes after each insert
  for Ins { pos, seq } in inss.iter().rev() {
    insert_at_inplace(&mut cons, *pos, &seq.chars().collect_vec());
  }

  cons.into_iter().filter(|&c| c != '\0').collect()
}

/// Updates the consensus sequence of the block and re-aligns the sequences to the new consensus.
fn update_block_consensus(block: &mut PangraphBlock, consensus: impl Into<String> + AsRef<str>) -> Result<(), Report> {
  // Reconstruct block sequences
  let seqs = block
    .alignments()
    .iter()
    .map(|(&nid, edit)| Ok((nid, edit.apply(block.consensus())?)))
    .collect::<Result<BTreeMap<NodeId, String>, Report>>()?;

  // debug assets: all sequences are non-empty
  #[cfg(any(debug_assertions, test))]
  {
    for (nid, seq) in &seqs {
      if seq.is_empty() {
        return make_error!(
          "node is empty!\nblock: {}\nnode: {}\nedits: {:?}\nconsensus: {}",
          block.id(),
          nid,
          block.alignments().get(nid),
          block.consensus()
        );
      }
    }
  }

  // Re-align sequences
  let alignments = seqs
    .into_iter()
    .map(|(nid, seq)| Ok((nid, map_variations(&consensus, &seq)?)))
    .collect::<Result<_, Report>>()?;

  *block = PangraphBlock::new(block.id(), consensus, alignments);

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::PangraphBlock;
  use crate::pangraph::pangraph_node::NodeId;
  use crate::pretty_assert_eq;
  use maplit::btreemap;

  fn block_1() -> PangraphBlock {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)    .T...--..........A.....
    //    node 2)    .T...--...C......|.....
    //    node 3)    .T...--...C.....--.....
    //    node 4)    .C.......---.....A.....
    //    node 5)    .....|...........A.....
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![Sub::new(1, 'T'), Sub::new(17, 'A')]),
      NodeId(2) => Edit::new(vec![Ins::new(17, "CAAT")],  vec![Del::new(5, 2)],     vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(3) => Edit::new(vec![],                      vec![Del::new(5, 2), Del::new(16,2)],     vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![],                      vec![Del::new(9, 3)],     vec![Sub::new(1, 'C'), Sub::new(17, 'A')]),
      NodeId(5) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![Sub::new(17, 'A')]),
    };
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_1_mut_reconsensus() -> PangraphBlock {
    let consensus = "ATGACTTCGATCTATTCAGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)    .....--................
    //    node 2)    .....--...C......G|....
    //    node 3)    .....--...C.....--.....
    //    node 4)    .C.......---...........
    //    node 5)    .G...|.................
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![]),
      NodeId(2) => Edit::new(vec![Ins::new(17, "CAAT")],  vec![Del::new(5, 2)],     vec![Sub::new(10, 'C'), Sub::new(17, 'G')]),
      NodeId(3) => Edit::new(vec![],                      vec![Del::new(5, 2), Del::new(16,2)],     vec![Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![],                      vec![Del::new(9, 3)],     vec![Sub::new(1, 'C')]),
      NodeId(5) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![Sub::new(1, 'G')]),
    };
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_2() -> PangraphBlock {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)   |.T.|.--......|...A..-..
    //    node 2)   |.T...--...C..|......--.|
    //    node 3)    .T..----..C............|
    //    node 4)    .C.|.....---.....A.....|
    //    node 5)   |...|.........|...A.--..
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![Ins::new(0, "G"), Ins::new(3, "AA"), Ins::new(13, "AA")],   vec![Del::new(5, 2), Del::new(20, 1)],  vec![Sub::new(1, 'T'), Sub::new(17, 'A')]),
      NodeId(2) => Edit::new(vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")],  vec![Del::new(5, 2), Del::new(20, 2)],  vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(3) => Edit::new(vec![Ins::new(23, "TT")],                                        vec![Del::new(4, 4)],                   vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![Ins::new(3, "C"), Ins::new(23, "TT")],                      vec![Del::new(9, 3)],                   vec![Sub::new(1, 'C'), Sub::new(17, 'A')]),
      NodeId(5) => Edit::new(vec![Ins::new(0, "G"), Ins::new(3, "C"), Ins::new(13, "AA")],    vec![Del::new(19, 2)],                  vec![Sub::new(17, 'A')])
    };
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_3() -> PangraphBlock {
    let consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC";
    let aln = btreemap! {
        NodeId(1) => Edit::new(vec![],                                        vec![Del::new(19, 4)],  vec![Sub::new(5, 'A')]),
        NodeId(2) => Edit::new(vec![Ins::new(35, "AA"), Ins::new(50, "TT")],  vec![Del::new(20, 3)],  vec![Sub::new(5, 'A')]),
        NodeId(3) => Edit::new(vec![],                                        vec![],                 vec![Sub::new(14, 'G'), Sub::new(27, 'G')]),
        NodeId(4) => Edit::new(vec![Ins::new(50, "TT")],                      vec![Del::new(20, 3)],  vec![Sub::new(5, 'A')]),
        NodeId(5) => Edit::new(vec![Ins::new(50, "TT")],                      vec![],                 vec![])
    };
    PangraphBlock::new(BlockId(3), consensus, aln)
  }

  fn block_3_reconsensus() -> PangraphBlock {
    let consensus = "GCCTCATCCCGACCACGCGTAACATGGGACAGGCCTGCGCTTGAGGCTT";
    let aln = btreemap! {
        NodeId(1) => Edit::new(vec![],                    vec![Del::new(19, 1), Del::new(47, 2)],   vec![]),
        NodeId(2) => Edit::new(vec![Ins::new(32, "AA")],  vec![],                                   vec![]),
        NodeId(3) => Edit::new(vec![Ins::new(20, "TAC")], vec![Del::new(47, 2)],                    vec![Sub::new(5, 'T'), Sub::new(14, 'G'), Sub::new(24, 'G')]),
        NodeId(4) => Edit::new(vec![],                    vec![],                                   vec![]),
        NodeId(5) => Edit::new(vec![Ins::new(20, "TAC")], vec![],                                   vec![Sub::new(5, 'T')])
    };
    PangraphBlock::new(BlockId(3), consensus, aln)
  }

  #[test]
  fn test_reconsensus_mutations() {
    let mut block = block_1();
    let expected_block = block_1_mut_reconsensus();
    reconsensus_mutations(&mut block).unwrap();
    pretty_assert_eq!(block, expected_block);
  }

  #[test]
  fn test_majority_deletions() {
    let dels = majority_deletions(&block_2());
    pretty_assert_eq!(dels, vec![5, 6, 20]);
  }

  #[test]
  fn test_majority_insertions() {
    let ins = majority_insertions(&block_2());
    pretty_assert_eq!(ins, vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")]);
  }

  #[test]
  fn test_apply_indels() {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    let dels = vec![5, 6, 20];
    let ins = vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")];
    let cons = apply_indels(consensus, &dels, &ins);
    pretty_assert_eq!(cons, "GAGGACCGATCTAAATTCGGAAATT");
  }

  #[test]
  fn test_reconsensus() {
    let mut block = block_3();
    let expected_block = block_3_reconsensus();
    reconsensus(&mut block).unwrap();
    pretty_assert_eq!(block.consensus(), expected_block.consensus());
    pretty_assert_eq!(block.alignments(), expected_block.alignments());
  }
}
