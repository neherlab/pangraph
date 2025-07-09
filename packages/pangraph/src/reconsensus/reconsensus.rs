use crate::align::map_variations::{BandParameters, map_variations};
use crate::commands::build::build_args::PangraphBuildArgs;
use crate::make_error;
use crate::pangraph::detach_unaligned::detach_unaligned_nodes;
use crate::pangraph::edits::{Del, Edit};
use crate::pangraph::edits::{Ins, Sub};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::NodeId;
use crate::utils::interval::positions_to_intervals;
// use crate::reconsensus::remove_nodes::remove_emtpy_nodes;
use crate::reconsensus::remove_nodes::find_empty_nodes;
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use rayon::prelude::*;
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
pub fn reconsensus_graph(
  graph: &mut Pangraph,
  ids_updated_blocks: Vec<BlockId>,
  args: &PangraphBuildArgs,
) -> Result<(), Report> {
  // there should be no empty nodes in the graph
  debug_assert!(
    find_empty_nodes(graph, &ids_updated_blocks).is_empty(),
    "Empty nodes found in the graph"
  );

  // reconsensus each block
  // i.e. update the consensus with majority variants
  // ids of blocks that undergo re-alignment are collected in realigned_block_ids
  let mut realigned_block_ids = Vec::new();
  for block_id in ids_updated_blocks {
    let block = graph.blocks.get_mut(&block_id).unwrap();
    let realigned = reconsensus(block, args)?;
    if realigned {
      realigned_block_ids.push(block_id);
    }
  }

  // For realigned blocks, pop them from graph.blocks, apply detach_unaligned_nodes to the list, and re-add them
  if !realigned_block_ids.is_empty() {
    let mut realigned_blocks = Vec::new();

    // Pop the realigned blocks from graph.blocks
    for block_id in &realigned_block_ids {
      if let Some(block) = graph.blocks.remove(block_id) {
        realigned_blocks.push(block);
      }
    }

    // Apply detach_unaligned_nodes. This removes unaligned nodes and re-adds them to the list as new blocks.
    detach_unaligned_nodes(&mut realigned_blocks, &mut graph.nodes)?;

    // Re-add all the blocks (including potentially new singleton blocks) to graph.blocks
    for block in realigned_blocks {
      graph.blocks.insert(block.id(), block);
    }
  }

  Ok(())
}

/// Performs the reconsensus operation inplace on a block.
/// Returns true or false, depending on whether sequences were re-aligned or not.
/// - if a position is mutated in > N/2 sites, it adds the mutation to the consensus and updates the alignment.
/// - if an in/del is present in > N/2 sites, it adds it to the consensus and re-aligns the sequences to the updated consensus.
fn reconsensus(block: &mut PangraphBlock, args: &PangraphBuildArgs) -> Result<bool, Report> {
  reconsensus_mutations(block)?;
  let ins = majority_insertions(block);
  let dels = majority_deletions(block);
  let re_align = !ins.is_empty() || !dels.is_empty();
  if re_align {
    let consensus = block.consensus();
    let consensus = apply_indels(consensus, &dels, &ins);
    // average shift of the new consensus with respect to the old one
    // once the indels are applied
    // used to help re-align the sequences
    let new_consensus_edits = Edit {
      subs: vec![],
      dels: positions_to_intervals(&dels)
        .iter()
        .map(|d| Del::new(d.start, d.len()))
        .collect(),
      inss: ins,
    };

    let new_consensus_band_params = BandParameters::from_edits(&new_consensus_edits, block.consensus_len())?;

    // debug assert: consensus is not empty
    debug_assert!(!consensus.is_empty(), "Consensus is empty after indels");

    update_block_consensus(block, &consensus, new_consensus_band_params, args)?;
  }
  Ok(re_align)
}

/// Re-computes the consensus for a block if a position is mutated in > N/2 sites.
fn reconsensus_mutations(block: &mut PangraphBlock) -> Result<(), Report> {
  let n = block.depth();
  let mut muts = btreemap! {};

  // count mutations
  for edit in block.alignments().values() {
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
    let original = block.consensus()[pos];

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
        },
        1 => {
          let s = &subs_at_pos[0];
          if s.alt == alt {
            edit.subs.retain(|sub| !(sub.pos == pos && sub.alt == alt));
          }
        },
        _ => {
          return make_error!(
            "At block {}: at position {pos}: sequence states disagree: {:}",
            block.id(),
            subs_at_pos
              .iter()
              .map(|sub| sub.alt.to_string())
              .collect_vec()
              .join(", ")
          );
        },
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
fn apply_indels(cons: impl Into<Seq>, dels: &[usize], inss: &[Ins]) -> Seq {
  let mut cons = cons.into();

  for &pos in dels {
    cons[pos] = AsciiChar(b'\0'); // Using '\0' to temporarily denote deleted positions
  }

  // Reverse to maintain correct insertion indexes after each insert
  for Ins { pos, seq } in inss.iter().rev() {
    cons.insert_seq(*pos, seq);
  }

  cons.retain(|&c| c != AsciiChar(b'\0'));

  cons
}

/// Updates the consensus sequence of the block and re-aligns the sequences to the new consensus.
fn update_block_consensus(
  block: &mut PangraphBlock,
  consensus: &Seq,
  new_consensus_band_params: BandParameters,
  args: &PangraphBuildArgs,
) -> Result<(), Report> {
  // Reconstruct block sequences
  let seqs = block
    .alignments()
    .iter()
    .map(|(&nid, edit)| {
      let seq = edit.apply(block.consensus())?;
      let old_band_params = BandParameters::from_edits(edit, block.consensus_len())?;
      let updated_band_params = BandParameters::new(
        old_band_params.mean_shift() - new_consensus_band_params.mean_shift(),
        old_band_params.band_width() + new_consensus_band_params.band_width(),
      );
      Ok((nid, (seq, updated_band_params)))
    })
    .collect::<Result<BTreeMap<NodeId, (Seq, BandParameters)>, Report>>()?;

  // debug assets: all sequences are non-empty
  #[cfg(any(debug_assertions, test))]
  {
    for (nid, (seq, _bp)) in &seqs {
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
    .into_par_iter()
    .map(|(nid, (seq, band_params))| Ok((nid, map_variations(consensus, &seq, band_params, args)?)))
    .collect::<Result<_, Report>>()?;

  *block = PangraphBlock::new(block.id(), consensus, alignments);

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
  use crate::pangraph::pangraph_node::NodeId;
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::pangraph_path::{PangraphPath, PathId};
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::utils::id::id;

  use maplit::btreemap;
  use pretty_assertions::assert_eq;

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

  fn block_mutations_only() -> PangraphBlock {
    let consensus = "ATGCGATCGATCGA";
    //               01234567890123
    //    node 1)    .C............  (mutation at pos 1: T->C)
    //    node 2)    .C............  (mutation at pos 1: T->C)
    //    node 3)    .C............  (mutation at pos 1: T->C)
    //    node 4)    ..........G...  (mutation at pos 10: C->G)
    //    node 5)    ..........G...  (mutation at pos 10: C->G)
    //     L = 14, N = 5
    // Position 1: T->C appears in 3/5 sequences (majority)
    // Position 10: C->G appears in 2/5 sequences (not majority)
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]),
      NodeId(2) => Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]),
      NodeId(3) => Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]),
      NodeId(4) => Edit::new(vec![], vec![], vec![Sub::new(10, 'G')]),
      NodeId(5) => Edit::new(vec![], vec![], vec![Sub::new(10, 'G')]),
    };
    PangraphBlock::new(BlockId(10), consensus, aln)
  }

  fn block_mutations_only_after_reconsensus() -> PangraphBlock {
    let consensus = "ACGCGATCGATCGA";
    //               01234567890123
    //    node 1)    ..............  (no mutations after consensus update)
    //    node 2)    ..............  (no mutations after consensus update)
    //    node 3)    ..............  (no mutations after consensus update)
    //    node 4)    .T........G...  (gets T at pos 1, keeps G at pos 10)
    //    node 5)    .T........G...  (gets T at pos 1, keeps G at pos 10)
    //     L = 14, N = 5
    // After reconsensus: consensus[1] = 'C' (majority), original T becomes mutation in nodes 4,5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![], vec![], vec![]),
      NodeId(2) => Edit::new(vec![], vec![], vec![]),
      NodeId(3) => Edit::new(vec![], vec![], vec![]),
      NodeId(4) => Edit::new(vec![], vec![], vec![Sub::new(1, 'T'), Sub::new(10, 'G')]),
      NodeId(5) => Edit::new(vec![], vec![], vec![Sub::new(1, 'T'), Sub::new(10, 'G')]),
    };
    PangraphBlock::new(BlockId(10), consensus, aln)
  }

  #[test]
  fn test_reconsensus_mutations() {
    let mut block = block_1();
    let expected_block = block_1_mut_reconsensus();
    reconsensus_mutations(&mut block).unwrap();
    assert_eq!(block, expected_block);
  }

  #[test]
  fn test_majority_deletions() {
    let dels = majority_deletions(&block_2());
    assert_eq!(dels, vec![5, 6, 20]);
  }

  #[test]
  fn test_majority_insertions() {
    let ins = majority_insertions(&block_2());
    assert_eq!(ins, vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")]);
  }

  #[test]
  fn test_apply_indels() {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    let dels = vec![5, 6, 20];
    let ins = vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")];
    let cons = apply_indels(consensus, &dels, &ins);
    assert_eq!(cons, "GAGGACCGATCTAAATTCGGAAATT");
  }

  #[test]
  fn test_reconsensus() {
    let mut block = block_3();
    let expected_block = block_3_reconsensus();
    let re_aligned = reconsensus(&mut block, &PangraphBuildArgs::default()).unwrap();
    assert!(re_aligned);
    assert_eq!(block.consensus(), expected_block.consensus());
    assert_eq!(block.alignments(), expected_block.alignments());
  }

  #[test]
  fn test_reconsensus_mutations_only_no_realignment() {
    let mut block = block_mutations_only();
    let expected_block = block_mutations_only_after_reconsensus();
    let re_aligned = reconsensus(&mut block, &PangraphBuildArgs::default()).unwrap();

    // Should return false because no indels require re-alignment
    assert!(!re_aligned);

    // But consensus should be updated and mutations healed
    assert_eq!(block.consensus(), expected_block.consensus());
    assert_eq!(block.alignments(), expected_block.alignments());
  }

  fn block_for_graph_test() -> PangraphBlock {
    let consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC";
    //               0         1         2         3         4
    //               01234567890123456789012345678901234567890123456789
    let aln = btreemap! {
        NodeId(1) => Edit::new(vec![],  vec![Del::new(0, 40)],   vec![]),  // deletion from position 0 to 40
        NodeId(2) => Edit::new(vec![],  vec![Del::new(35, 15)],  vec![]),  // deletion from position 35 to end (50-35=15)
        NodeId(3) => Edit::new(vec![],  vec![Del::new(35, 15)],  vec![]),  // deletion from position 35 to end
        NodeId(4) => Edit::new(vec![],  vec![Del::new(35, 15)],  vec![]),  // deletion from position 35 to end
        NodeId(5) => Edit::new(vec![],  vec![],                  vec![])   // no deletions
    };
    PangraphBlock::new(BlockId(20), consensus, aln)
  }

  fn block_for_graph_test_expected() -> PangraphBlock {
    // After reconsensus, the majority deletions (positions 35-49) should be applied
    // Original consensus: "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC" (50 chars)
    // Expected consensus: "GCCTCTTCCCGACCACGCGTTACAACATGGGACAG" (35 chars)
    let consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAG";
    let aln = btreemap! {
        NodeId(2) => Edit::new(vec![],                           vec![],                  vec![]),       // no edits (was majority deletion)
        NodeId(3) => Edit::new(vec![],                           vec![],                  vec![]),       // no edits (was majority deletion)
        NodeId(4) => Edit::new(vec![],                           vec![],                  vec![]),       // no edits (was majority deletion)
        NodeId(5) => Edit::new(vec![Ins::new(35, "GCCTGCGCTTGAGGC")], vec![],                  vec![])        // insertion to represent chars that were deleted from consensus
    };
    // NodeId(1) was detached because it was empty after deletions
    PangraphBlock::new(BlockId(20), consensus, aln)
  }

  fn signleton_block_expected() -> PangraphBlock {
    // This is the singleton block that should be created for NodeId(1) after reconsensus
    // this should be the reverse-complement of CGCTTGAGGC, because NodeId(1) was in Reverse strand
    let consensus = Seq::from("GCCTCAAGCG");
    PangraphBlock::from_consensus(consensus.clone(), BlockId(id((NodeId(1), &consensus))), NodeId(1))
  }

  #[test]
  fn test_reconsensus_graph() {
    // Create a block that will be modified by reconsensus
    let initial_block = block_for_graph_test();
    let expected_block = block_for_graph_test_expected();
    let singleton_block_exp = signleton_block_expected();

    // Create nodes for the block with lengths reflecting actual sequence lengths
    let nodes = btreemap! {
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), initial_block.id(), PathId(1), Reverse, (0, 10)),   // 50 - 40 = 9 (deletes positions 0-39)
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), initial_block.id(), PathId(2), Forward, (0, 35)),  // 50 - 15 = 35 (deletes positions 35-49)
      NodeId(3) => PangraphNode::new(Some(NodeId(3)), initial_block.id(), PathId(3), Forward, (0, 35)),  // 50 - 15 = 35 (deletes positions 35-49)
      NodeId(4) => PangraphNode::new(Some(NodeId(4)), initial_block.id(), PathId(4), Forward, (0, 35)),  // 50 - 15 = 35 (deletes positions 35-49)
      NodeId(5) => PangraphNode::new(Some(NodeId(5)), initial_block.id(), PathId(5), Forward, (0, 49)),  // no deletions
    };

    // Create paths
    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1)], 49, false, None, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [NodeId(2)], 49, false, None, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [NodeId(3)], 49, false, None, None),
      PathId(4) => PangraphPath::new(Some(PathId(4)), [NodeId(4)], 49, false, None, None),
      PathId(5) => PangraphPath::new(Some(PathId(5)), [NodeId(5)], 49, false, None, None),
    };

    // Create blocks map
    let blocks = btreemap! {
      initial_block.id() => initial_block.clone()
    };

    // Create the graph
    let mut graph = Pangraph { paths, blocks, nodes };

    // Apply reconsensus_graph
    let result = reconsensus_graph(&mut graph, vec![initial_block.id()], &PangraphBuildArgs::default());

    // Check that the operation succeeded
    result.unwrap();

    // Get the two created blocks, the modified and the singleton
    let final_block = &graph.blocks[&initial_block.id()];
    let singleton_block = &graph.blocks[&singleton_block_exp.id()];

    // Direct comparison with expected result
    assert_eq!(final_block.consensus(), expected_block.consensus());
    assert_eq!(final_block.alignments(), expected_block.alignments());

    // Check that the empty node was detached and a new singleton block was created
    assert_eq!(singleton_block.consensus(), singleton_block_exp.consensus());
    assert_eq!(singleton_block.alignments(), singleton_block_exp.alignments());

    // check that the node was updated correctly, flipping the strandedness
    let new_node1 = &graph.nodes[&NodeId(1)];
    let expected_node1 = PangraphNode::new(Some(NodeId(1)), singleton_block_exp.id(), PathId(1), Forward, (0, 10));
    assert_eq!(new_node1, &expected_node1);
  }
}
