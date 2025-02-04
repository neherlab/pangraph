#![allow(non_snake_case)]

use crate::pangraph::edits::{Del, Edit, Ins, Sub};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::PangraphBlock;
use crate::pangraph::pangraph_interval::PangraphInterval;
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::pangraph::strand::Strand;
use crate::representation::seq::Seq;
use std::collections::BTreeMap;

pub fn slice_substitutions(i: &PangraphInterval, S: &[Sub]) -> Vec<Sub> {
  S.iter()
    .filter(|s| i.contains(s.pos))
    .map(|s| Sub {
      pos: s.pos - i.interval.start,
      alt: s.alt,
    })
    .collect()
}

pub fn slice_deletions(i: &PangraphInterval, D: &[Del]) -> Vec<Del> {
  D.iter()
    .filter(|d| i.has_overlap_with(&d.interval()))
    .map(|d| {
      let new_start = d.pos.max(i.interval.start) - i.interval.start;
      let new_end = (d.pos + d.len).min(i.interval.end) - i.interval.start;
      let new_len = new_end - new_start;
      Del {
        pos: new_start,
        len: new_len,
      }
    })
    .collect()
}

pub fn slice_insertions(i: &PangraphInterval, I: &[Ins], block_L: usize) -> Vec<Ins> {
  I.iter()
    .filter(|ins| i.insertion_overlap(ins.pos, block_L))
    .map(|ins| Ins {
      pos: ins.pos - i.interval.start,
      seq: ins.seq.clone(),
    })
    .collect()
}

pub fn slice_edits(i: &PangraphInterval, ed: &Edit, block_L: usize) -> Edit {
  Edit {
    inss: slice_insertions(i, &ed.inss, block_L),
    dels: slice_deletions(i, &ed.dels),
    subs: slice_substitutions(i, &ed.subs),
  }
}

pub fn new_strandedness(old_strandedness: Strand, orientation: Strand, is_anchor: bool) -> Strand {
  if is_anchor || orientation.is_forward() {
    old_strandedness
  } else {
    old_strandedness.reverse()
  }
}

// given the old position of a node, the coordinates of the node in the new block,
// the length of the path and the strandedness of the node, return the new node position.
// Accounts for circular paths.
// Nb: a node that spans the whole path [0, L] will have position [0, 0]
pub fn new_position_circular(
  old_position: (usize, usize),
  node_coords: (usize, usize),
  path_L: usize,
  old_strandedness: Strand,
) -> (usize, usize) {
  let (old_s, old_e) = old_position;
  let (new_s_in_node, new_e_in_node) = node_coords;
  if old_strandedness.is_forward() {
    ((old_s + new_s_in_node) % path_L, (old_s + new_e_in_node) % path_L)
  } else {
    (
      (old_e + path_L - new_e_in_node) % path_L,
      (old_e + path_L - new_s_in_node) % path_L,
    )
  }
}

// given the old position of a node, the coordinates of the node in the new block,
// and the strandedness of the node, return the new node position.
// Does not account for circular paths.
// Nb: a node that spans the whole path [0, L] will have position [0, L]
pub fn new_position_non_circular(
  old_position: (usize, usize),
  node_coords: (usize, usize),
  old_strandedness: Strand,
) -> (usize, usize) {
  let (old_s, old_e) = old_position;
  let (new_s_in_node, new_e_in_node) = node_coords;
  if old_strandedness.is_forward() {
    (old_s + new_s_in_node, old_s + new_e_in_node)
  } else {
    (old_e - new_e_in_node, old_e - new_s_in_node)
  }
}

pub fn interval_node_coords(i: &PangraphInterval, edits: &Edit, block_L: usize) -> (usize, usize) {
  let (mut s, mut e) = (i.interval.start, i.interval.end);
  for d in &edits.dels {
    if d.pos <= i.interval.start {
      let l = (d.len + d.pos).min(i.interval.start) - d.pos;
      s -= l;
    }
    if d.pos < i.interval.end {
      let l = (d.len + d.pos).min(i.interval.end) - d.pos;
      e -= l;
    }
  }
  for ins in &edits.inss {
    if ins.pos < i.interval.start {
      s += ins.seq.len();
    }
    if ins.pos < i.interval.end {
      e += ins.seq.len();
    }
    if ins.pos == i.interval.end && ins.pos == block_L {
      e += ins.seq.len();
    }
  }
  (s, e)
}

/// given a block, an interval and the graph, it extract the slice of the block
/// defined by the interval.
///
/// It returns the new block, along with a dictionary of updates for the nodes,
/// that maps old node ids to new nodes.
/// In case a slice contains an empty node, the corresponding entry in the dictionary
/// will be None.
pub fn block_slice(
  b: &PangraphBlock,
  i: &PangraphInterval,
  G: &Pangraph,
) -> (PangraphBlock, BTreeMap<NodeId, Option<PangraphNode>>) {
  #[allow(clippy::string_slice)]
  // consensus of the new block
  let new_consensus: Seq = b.consensus()[i.interval.to_range()].into();

  // length of the new block
  let block_L = b.consensus_len();
  debug_assert!(block_L > 0, "Block {} has length 0", b.id());

  // containers for new alignment and node updates
  let mut node_updates = BTreeMap::new();
  let mut new_alignment = BTreeMap::new();

  for (old_node_id, edits) in b.alignments() {
    #[cfg(any(debug_assertions, test))]
    edits.sanity_check(b.consensus_len()).unwrap();

    let old_node = &G.nodes[old_node_id];
    let old_strandedness = old_node.strand();

    let new_strand = if i.aligned {
      new_strandedness(
        old_strandedness,
        i.orientation.unwrap(), // FIXME
        i.is_anchor.unwrap(),   // FIXME
      )
    } else {
      old_strandedness
    };

    let path_L = G.paths[&old_node.path_id()].tot_len;
    let node_coords = interval_node_coords(i, edits, block_L);
    let circular = G.paths[&old_node.path_id()].circular();
    let new_pos = if circular {
      new_position_circular(old_node.position(), node_coords, path_L, old_strandedness)
    } else {
      new_position_non_circular(old_node.position(), node_coords, old_strandedness)
    };

    let new_node = PangraphNode::new(None, i.new_block_id, old_node.path_id(), new_strand, new_pos);

    // extract edits for the slice
    let new_edits = slice_edits(i, edits, block_L);

    #[cfg(any(debug_assertions, test))]
    new_edits.sanity_check(new_consensus.len()).unwrap();

    if new_edits.is_empty_alignment(&new_consensus) {
      // if the node is empty, add `None` to the node updates
      // and do not add the alignment to the new block
      node_updates.insert(*old_node_id, None);
    } else {
      // if the node is not empty, add the alignment to the new block
      let ovw = new_alignment.insert(new_node.id(), new_edits);
      debug_assert!(ovw.is_none(), "Node id was already present! {ovw:?}");
      node_updates.insert(*old_node_id, Some(new_node));
    }
  }

  let new_block = PangraphBlock::new(i.new_block_id, new_consensus, new_alignment);

  (new_block, node_updates)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::BlockId;
  use crate::pangraph::pangraph_interval::*;
  use crate::pangraph::pangraph_node::*;
  use crate::pangraph::pangraph_path::{PangraphPath, PathId};
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::utils::interval::Interval;
  use maplit::btreemap;

  fn generate_example() -> (String, Edit) {
    let seq = o!("ACTGGATATCCGATATTCGAG");
    let ed = Edit {
      subs: vec![
        Sub::new(2, 'C'),
        Sub::new(5, 'C'),
        Sub::new(6, 'G'),
        Sub::new(7, 'C'),
        Sub::new(13, 'G'),
        Sub::new(14, 'T'),
        Sub::new(18, 'C'),
        Sub::new(20, 'A'),
      ],
      dels: vec![
        Del::new(0, 2),
        Del::new(4, 3),
        Del::new(9, 2),
        Del::new(13, 4),
        Del::new(18, 3),
      ],
      inss: vec![
        Ins::new(2, "CC"),
        Ins::new(5, "A"),
        Ins::new(6, "TTT"),
        Ins::new(10, "C"),
        Ins::new(13, "T"),
        Ins::new(14, "GG"),
        Ins::new(17, "A"),
        Ins::new(21, "A"),
      ],
    };
    (seq, ed)
  }

  #[test]
  fn test_slice_substitutions() {
    let (_, ed) = generate_example();
    let i = PangraphInterval {
      interval: Interval { start: 6, end: 14 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_substitutions(&i, &ed.subs);
    assert_eq!(new_ed, vec![Sub::new(0, 'G'), Sub::new(1, 'C'), Sub::new(7, 'G'),]);
    let i = PangraphInterval {
      interval: Interval { start: 15, end: 21 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_substitutions(&i, &ed.subs);
    assert_eq!(new_ed, vec![Sub::new(3, 'C'), Sub::new(5, 'A')]);
  }

  #[test]
  fn test_slice_deletions() {
    let (_, ed) = generate_example();
    let i = PangraphInterval {
      interval: Interval { start: 6, end: 14 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_deletions(&i, &ed.dels);
    assert_eq!(
      new_ed,
      vec![Del { pos: 0, len: 1 }, Del { pos: 3, len: 2 }, Del { pos: 7, len: 1 },]
    );
    let i = PangraphInterval {
      interval: Interval { start: 15, end: 21 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_deletions(&i, &ed.dels);
    assert_eq!(new_ed, vec![Del { pos: 0, len: 2 }, Del { pos: 3, len: 3 },]);
  }

  #[test]
  fn test_slice_insertions() {
    let (seq, ed) = generate_example();
    let i = PangraphInterval {
      interval: Interval { start: 6, end: 14 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_insertions(&i, &ed.inss, seq.len());
    assert_eq!(new_ed, vec![Ins::new(0, "TTT"), Ins::new(4, "C"), Ins::new(7, "T"),]);
    let i = PangraphInterval {
      interval: Interval { start: 15, end: 21 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_insertions(&i, &ed.inss, seq.len());
    assert_eq!(new_ed, vec![Ins::new(2, "A"), Ins::new(6, "A")]);
  }

  #[test]
  fn test_interval_node_coords() {
    let (seq, ed) = generate_example();

    let i = PangraphInterval {
      interval: Interval { start: 6, end: 14 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_pos = interval_node_coords(&i, &ed, seq.len());
    assert_eq!(new_pos, (5, 14));

    let i = PangraphInterval {
      interval: Interval { start: 15, end: 21 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_pos = interval_node_coords(&i, &ed, seq.len());
    assert_eq!(new_pos, (16, 19));
  }

  #[test]
  fn test_new_position_circular() {
    let path_L = 100;

    let strandedness = Forward;
    let node_coords = (10, 20);
    let old_position = (10, 40);
    let new_pos = new_position_circular(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (20, 30));

    let old_position = (95, 20);
    let new_pos = new_position_circular(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (5, 15));

    let strandedness = Reverse;
    let old_position = (10, 50);
    let new_pos = new_position_circular(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (30, 40));

    let old_position = (40, 5);
    let new_pos = new_position_circular(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (85, 95));

    let strandedness = Forward;
    let old_position = (0, 100);
    let node_coords = (0, 100);
    let new_pos = new_position_circular(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (0, 0));
  }

  #[test]
  fn test_new_position_non_circular() {
    let strandedness = Forward;
    let node_coords = (10, 20);
    let old_position = (10, 40);
    let new_pos = new_position_non_circular(old_position, node_coords, strandedness);
    assert_eq!(new_pos, (20, 30));

    let strandedness = Reverse;
    let node_coords = (10, 20);
    let old_position = (10, 50);
    let new_pos = new_position_non_circular(old_position, node_coords, strandedness);
    assert_eq!(new_pos, (30, 40));

    let strandedness = Forward;
    let node_coords = (0, 10);
    let old_position = (0, 20);
    let new_pos = new_position_non_circular(old_position, node_coords, strandedness);
    assert_eq!(new_pos, (0, 10));

    let strandedness = Forward;
    let node_coords = (0, 100);
    let old_position = (0, 100);
    let new_pos = new_position_non_circular(old_position, node_coords, strandedness);
    assert_eq!(new_pos, (0, 100));
  }

  #[test]
  fn test_node_coords() {
    let i = PangraphInterval {
      interval: Interval { start: 10, end: 20 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let ed = Edit {
      subs: vec![Sub::new(2, 'G'), Sub::new(13, 'T'), Sub::new(24, 'T')],
      dels: vec![Del { pos: 18, len: 3 }],
      inss: vec![Ins::new(7, "A"), Ins::new(10, "AAAA"), Ins::new(20, "TTTTTTTT")],
    };
    let block_L = 100;
    let new_pos = interval_node_coords(&i, &ed, block_L);
    assert_eq!(new_pos, (11, 23));
  }

  #[test]
  fn test_block_slice_fwd_anchor() {
    let (b, G) = generate_block_example();
    let new_bid = BlockId(42);

    let i = PangraphInterval {
      interval: Interval::new(10, 20),
      aligned: true,
      new_block_id: new_bid,
      is_anchor: Some(true),
      orientation: Some(Forward),
    };

    let (new_b, new_nodes) = block_slice(&b, &i, &G);
    assert_eq!(new_b.consensus(), "TATATTTATC");

    let nn1 = PangraphNode::new(None, new_bid, PathId(1), Forward, (111, 120));
    let nn1_slice = new_nodes[&NodeId(1)].as_ref().unwrap();
    assert_eq!(nn1, *nn1_slice);

    let nn2 = PangraphNode::new(None, new_bid, PathId(2), Reverse, (1008, 1017));
    let nn2_slice = new_nodes[&NodeId(2)].as_ref().unwrap();
    assert_eq!(nn2, *nn2_slice);

    let nn3 = PangraphNode::new(None, new_bid, PathId(3), Reverse, (96, 4));
    let nn3_slice = new_nodes[&NodeId(3)].as_ref().unwrap();
    assert_eq!(nn3, *nn3_slice);

    assert_eq!(
      new_nodes,
      btreemap! {
        NodeId(1) => Some(nn1.clone()),
        NodeId(2) => Some(nn2.clone()),
        NodeId(3) => Some(nn3.clone()),
      }
    );

    let n1ed = Edit {
      subs: vec![Sub::new(3, 'T')],
      dels: vec![Del::new(8, 2)],
      inss: vec![Ins::new(0, "A")],
    };
    assert_eq!(new_b.alignment(nn1.id()), &n1ed);

    let n2ed = Edit {
      subs: vec![Sub::new(9, 'G')],
      dels: vec![Del::new(3, 2)],
      inss: vec![Ins::new(7, "T")],
    };
    assert_eq!(new_b.alignment(nn2.id()), &n2ed);

    let n3ed = Edit {
      subs: vec![],
      dels: vec![Del { pos: 0, len: 2 }],
      inss: vec![],
    };
    assert_eq!(new_b.alignment(nn3.id()), &n3ed);
  }

  fn generate_block_example() -> (PangraphBlock, Pangraph) {
    let seq = o!("ACTTGATCCTTATATTTATCCGATCAT");
    let bid = BlockId(1);

    let ed1 = Edit {
      subs: vec![Sub::new(2, 'G'), Sub::new(13, 'T'), Sub::new(24, 'T')],
      dels: vec![Del { pos: 18, len: 3 }],
      inss: vec![Ins::new(7, "A"), Ins::new(10, "A")],
    };

    let ed2 = Edit {
      subs: vec![Sub::new(4, 'T'), Sub::new(19, 'G'), Sub::new(20, 'G')],
      dels: vec![Del { pos: 6, len: 2 }, Del { pos: 13, len: 2 }],
      inss: vec![Ins::new(17, "T"), Ins::new(25, "A")],
    };

    let ed3 = Edit {
      subs: vec![],
      dels: vec![Del { pos: 2, len: 4 }, Del { pos: 9, len: 3 }, Del { pos: 24, len: 2 }],
      inss: vec![Ins::new(20, "T")],
    };

    let n1 = PangraphNode::new(Some(NodeId(1)), bid, PathId(1), Forward, (100, 125));
    let n2 = PangraphNode::new(Some(NodeId(2)), bid, PathId(2), Reverse, (1000, 1025));
    let n3 = PangraphNode::new(Some(NodeId(3)), bid, PathId(3), Reverse, (90, 9));

    let p1 = PangraphPath::new(Some(PathId(1)), /*"p1"*/ [NodeId(1), NodeId(4)], 2000, true, None, None);
    let p2 = PangraphPath::new(Some(PathId(2)), /*"p2"*/ [NodeId(2), NodeId(5)], 2000, true, None, None);
    let p3 = PangraphPath::new(Some(PathId(3)), /*"p3"*/ [NodeId(3), NodeId(6)], 100, true, None, None);

    let b1 = PangraphBlock::new(
      bid,
      Seq::from_str(&seq),
      btreemap! {
        NodeId(1) => ed1,
        NodeId(2) => ed2,
        NodeId(3) => ed3,
      },
    );

    let G = Pangraph {
      paths: btreemap! {
        PathId(1) => p1,
        PathId(2) => p2,
        PathId(3) => p3,
      },
      blocks: btreemap! {
        bid => b1.clone(),
      },
      nodes: btreemap! {
        NodeId(1) => n1,
        NodeId(2) => n2,
        NodeId(3) => n3,
      },
    };

    (b1, G)
  }

  #[test]
  fn test_block_slice_rev_append() {
    let (b, G) = generate_block_example();
    let new_bid = BlockId(42);
    let i = PangraphInterval {
      interval: Interval { start: 10, end: 20 },
      aligned: true,
      new_block_id: new_bid,
      is_anchor: Some(false),
      orientation: Some(Reverse),
    };

    let (new_b, new_nodes) = block_slice(&b, &i, &G);

    assert_eq!(new_b.consensus(), "TATATTTATC");

    let nn1 = PangraphNode::new(None, new_bid, PathId(1), Reverse, (111, 120));
    let nn1_slice = new_nodes[&NodeId(1)].as_ref().unwrap();
    assert_eq!(nn1, *nn1_slice);

    let nn2 = PangraphNode::new(None, new_bid, PathId(2), Forward, (1008, 1017));
    let nn2_slice = new_nodes[&NodeId(2)].as_ref().unwrap();
    assert_eq!(nn2, *nn2_slice);

    let nn3 = PangraphNode::new(None, new_bid, PathId(3), Forward, (96, 4));
    let nn3_slice = new_nodes[&NodeId(3)].as_ref().unwrap();
    assert_eq!(nn3, *nn3_slice);

    assert_eq!(
      new_nodes,
      btreemap! {
        NodeId(1) => Some(nn1.clone()),
        NodeId(2) => Some(nn2.clone()),
        NodeId(3) => Some(nn3.clone()),
      }
    );

    let n1ed = Edit {
      subs: vec![Sub::new(3, 'T')],
      dels: vec![Del { pos: 8, len: 2 }],
      inss: vec![Ins::new(0, "A")],
    };
    assert_eq!(new_b.alignment(nn1.id()), &n1ed);

    let n2ed = Edit {
      subs: vec![Sub::new(9, 'G')],
      dels: vec![Del { pos: 3, len: 2 }],
      inss: vec![Ins::new(7, "T")],
    };
    assert_eq!(new_b.alignment(nn2.id()), &n2ed);

    let n3ed = Edit {
      subs: vec![],
      dels: vec![Del { pos: 0, len: 2 }],
      inss: vec![],
    };
    assert_eq!(new_b.alignment(nn3.id()), &n3ed);
  }
}
