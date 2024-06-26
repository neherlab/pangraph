#![allow(non_snake_case)]

use crate::pangraph::edits::{Del, Edit, Ins, Sub};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::PangraphBlock;
use crate::pangraph::pangraph_interval::PangraphInterval;
use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
use crate::utils::id::Id;
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

pub fn new_strandedness(old_strandedness: bool, orientation: bool, is_anchor: bool) -> bool {
  if is_anchor || orientation {
    old_strandedness
  } else {
    !old_strandedness
  }
}

pub fn new_position(
  old_position: (usize, usize),
  node_coords: (usize, usize),
  path_L: usize,
  old_strandedness: bool,
) -> (usize, usize) {
  let (old_s, old_e) = old_position;
  let (new_s_in_node, new_e_in_node) = node_coords;
  if old_strandedness {
    ((old_s + new_s_in_node) % path_L, (old_s + new_e_in_node) % path_L)
  } else {
    (
      (old_e + path_L - new_e_in_node) % path_L,
      (old_e + path_L - new_s_in_node) % path_L,
    )
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

pub fn block_slice(
  b: &PangraphBlock,
  i: &PangraphInterval,
  G: &Pangraph,
) -> (PangraphBlock, BTreeMap<NodeId, PangraphNode>) {
  #[allow(clippy::string_slice)]
  let new_consensus = b.consensus[i.interval.to_range()].to_owned();
  let block_L = b.consensus_len();

  let mut node_updates = BTreeMap::new();
  let mut new_alignment = BTreeMap::new();

  for (old_node_id, edits) in &b.alignments {
    let old_node = &G.nodes[old_node_id];
    let old_strandedness = old_node.strand();

    let new_strand = if i.aligned {
      new_strandedness(
        old_strandedness,
        i.orientation.unwrap_or(false),
        i.is_anchor.unwrap_or(false),
      )
    } else {
      old_strandedness
    };

    let path_L = G.paths[&old_node.path_id()].tot_len;
    let node_coords = interval_node_coords(i, edits, block_L);
    let new_pos = new_position(old_node.position(), node_coords, path_L, old_strandedness);

    let new_node = PangraphNode::new(i.new_block_id, old_node.path_id(), new_strand, new_pos);
    node_updates.insert(*old_node_id, new_node.clone());

    let new_edits = slice_edits(i, edits, block_L);
    new_alignment.insert(new_node.id(), new_edits);
  }

  let new_block = PangraphBlock {
    consensus: new_consensus,
    alignments: new_alignment,
  };

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
  use crate::utils::interval::Interval;
  use maplit::btreemap;

  fn generate_example() -> (String, Edit) {
    let seq = o!("ACTGGATATCCGATATTCGAG");
    let ed = Edit {
      subs: vec![
        Sub { pos: 2, alt: 'C' },
        Sub { pos: 5, alt: 'C' },
        Sub { pos: 6, alt: 'G' },
        Sub { pos: 7, alt: 'C' },
        Sub { pos: 13, alt: 'G' },
        Sub { pos: 14, alt: 'T' },
        Sub { pos: 18, alt: 'C' },
        Sub { pos: 20, alt: 'A' },
      ],
      dels: vec![
        Del { pos: 0, len: 2 },
        Del { pos: 4, len: 3 },
        Del { pos: 9, len: 2 },
        Del { pos: 13, len: 4 },
        Del { pos: 18, len: 3 },
      ],
      inss: vec![
        Ins { pos: 2, seq: o!("CC") },
        Ins { pos: 5, seq: o!("A") },
        Ins { pos: 6, seq: o!("TTT") },
        Ins { pos: 10, seq: o!("C") },
        Ins { pos: 13, seq: o!("T") },
        Ins { pos: 14, seq: o!("GG") },
        Ins { pos: 17, seq: o!("A") },
        Ins { pos: 21, seq: o!("A") },
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
    assert_eq!(
      new_ed,
      vec![
        Sub { pos: 0, alt: 'G' },
        Sub { pos: 1, alt: 'C' },
        Sub { pos: 7, alt: 'G' },
      ]
    );
    let i = PangraphInterval {
      interval: Interval { start: 15, end: 21 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_substitutions(&i, &ed.subs);
    assert_eq!(new_ed, vec![Sub { pos: 3, alt: 'C' }, Sub { pos: 5, alt: 'A' },]);
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
    assert_eq!(
      new_ed,
      vec![
        Ins { pos: 0, seq: o!("TTT") },
        Ins { pos: 4, seq: o!("C") },
        Ins { pos: 7, seq: o!("T") },
      ]
    );
    let i = PangraphInterval {
      interval: Interval { start: 15, end: 21 },
      aligned: true,
      new_block_id: BlockId(0),
      is_anchor: None,
      orientation: None,
    };
    let new_ed = slice_insertions(&i, &ed.inss, seq.len());
    assert_eq!(
      new_ed,
      vec![Ins { pos: 2, seq: o!("A") }, Ins { pos: 6, seq: o!("A") },]
    );
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
  fn test_new_position() {
    let path_L = 100;

    let strandedness = true;
    let node_coords = (10, 20);
    let old_position = (10, 40);
    let new_pos = new_position(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (20, 30));

    let old_position = (95, 20);
    let new_pos = new_position(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (5, 15));

    let strandedness = false;
    let old_position = (10, 50);
    let new_pos = new_position(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (30, 40));

    let old_position = (40, 5);
    let new_pos = new_position(old_position, node_coords, path_L, strandedness);
    assert_eq!(new_pos, (85, 95));
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
      subs: vec![
        Sub { pos: 2, alt: 'G' },
        Sub { pos: 13, alt: 'T' },
        Sub { pos: 24, alt: 'T' },
      ],
      dels: vec![Del { pos: 18, len: 3 }],
      inss: vec![
        Ins { pos: 7, seq: o!("A") },
        Ins {
          pos: 10,
          seq: o!("AAAA"),
        },
        Ins {
          pos: 20,
          seq: o!("TTTTTTTT"),
        },
      ],
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
      orientation: Some(true),
    };

    let (new_b, new_nodes) = block_slice(&b, &i, &G);
    assert_eq!(new_b.consensus, "TATATTTATC");

    let nn1 = PangraphNode::new(new_bid, PathId(1), true, (111, 120));
    assert_eq!(&new_nodes[&NodeId(1)], &nn1);

    let nn2 = PangraphNode::new(new_bid, PathId(2), false, (1008, 1017));
    assert_eq!(&new_nodes[&NodeId(2)], &nn2);

    let nn3 = PangraphNode::new(new_bid, PathId(3), false, (96, 4));
    assert_eq!(&new_nodes[&NodeId(3)], &nn3);

    assert_eq!(
      new_nodes,
      btreemap! {
        NodeId(1) => nn1.clone(),
        NodeId(2) => nn2.clone(),
        NodeId(3) => nn3.clone(),
      }
    );

    let n1ed = Edit {
      subs: vec![Sub { pos: 3, alt: 'T' }],
      dels: vec![Del { pos: 8, len: 2 }],
      inss: vec![Ins { pos: 0, seq: o!("A") }],
    };
    assert_eq!(&new_b.alignments[&nn1.id()], &n1ed);

    let n2ed = Edit {
      subs: vec![Sub { pos: 9, alt: 'G' }],
      dels: vec![Del { pos: 3, len: 2 }],
      inss: vec![Ins { pos: 7, seq: o!("T") }],
    };
    assert_eq!(&new_b.alignments[&nn2.id()], &n2ed);

    let n3ed = Edit {
      subs: vec![],
      dels: vec![Del { pos: 0, len: 2 }],
      inss: vec![],
    };
    assert_eq!(&new_b.alignments[&nn3.id()], &n3ed);
  }

  fn generate_block_example() -> (PangraphBlock, Pangraph) {
    let seq = o!("ACTTGATCCTTATATTTATCCGATCAT");
    let bid = BlockId(1);

    let ed1 = Edit {
      subs: vec![
        Sub { pos: 2, alt: 'G' },
        Sub { pos: 13, alt: 'T' },
        Sub { pos: 24, alt: 'T' },
      ],
      dels: vec![Del { pos: 18, len: 3 }],
      inss: vec![Ins { pos: 7, seq: o!("A") }, Ins { pos: 10, seq: o!("A") }],
    };

    let ed2 = Edit {
      subs: vec![
        Sub { pos: 4, alt: 'T' },
        Sub { pos: 19, alt: 'G' },
        Sub { pos: 20, alt: 'G' },
      ],
      dels: vec![Del { pos: 6, len: 2 }, Del { pos: 13, len: 2 }],
      inss: vec![Ins { pos: 17, seq: o!("T") }, Ins { pos: 25, seq: o!("A") }],
    };

    let ed3 = Edit {
      subs: vec![],
      dels: vec![Del { pos: 2, len: 4 }, Del { pos: 9, len: 3 }, Del { pos: 24, len: 2 }],
      inss: vec![Ins { pos: 20, seq: o!("T") }],
    };

    let n1 = PangraphNode::new(bid, PathId(1), true, (100, 125));
    let n2 = PangraphNode::new(bid, PathId(2), false, (1000, 1025));
    let n3 = PangraphNode::new(bid, PathId(3), false, (90, 9));

    let mut p1 = PangraphPath::new("p1", &[NodeId(1), NodeId(4)], true);
    p1.tot_len = 2000;

    let mut p2 = PangraphPath::new("p2", &[NodeId(2), NodeId(5)], true);
    p2.tot_len = 2000;

    let mut p3 = PangraphPath::new("p3", &[NodeId(3), NodeId(6)], true);
    p3.tot_len = 100;

    let b1 = PangraphBlock {
      consensus: seq,
      alignments: btreemap! {
        NodeId(1) => ed1,
        NodeId(2) => ed2,
        NodeId(3) => ed3,
      },
    };

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
      orientation: Some(false),
    };

    let (new_b, new_nodes) = block_slice(&b, &i, &G);

    assert_eq!(new_b.consensus, "TATATTTATC");

    let nn1 = PangraphNode::new(new_bid, PathId(1), false, (111, 120));
    assert_eq!(&new_nodes[&NodeId(1)], &nn1);

    let nn2 = PangraphNode::new(new_bid, PathId(2), true, (1008, 1017));
    assert_eq!(&new_nodes[&NodeId(2)], &nn2);

    let nn3 = PangraphNode::new(new_bid, PathId(3), true, (96, 4));
    assert_eq!(&new_nodes[&NodeId(3)], &nn3);

    assert_eq!(
      new_nodes,
      btreemap! {
        NodeId(1) => nn1.clone(),
        NodeId(2) => nn2.clone(),
        NodeId(3) => nn3.clone(),
      }
    );

    let n1ed = Edit {
      subs: vec![Sub { pos: 3, alt: 'T' }],
      dels: vec![Del { pos: 8, len: 2 }],
      inss: vec![Ins { pos: 0, seq: o!("A") }],
    };
    assert_eq!(&new_b.alignments[&nn1.id()], &n1ed);

    let n2ed = Edit {
      subs: vec![Sub { pos: 9, alt: 'G' }],
      dels: vec![Del { pos: 3, len: 2 }],
      inss: vec![Ins { pos: 7, seq: o!("T") }],
    };
    assert_eq!(&new_b.alignments[&nn2.id()], &n2ed);

    let n3ed = Edit {
      subs: vec![],
      dels: vec![Del { pos: 0, len: 2 }],
      inss: vec![],
    };
    assert_eq!(&new_b.alignments[&nn3.id()], &n3ed);
  }
}
