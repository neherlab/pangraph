use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_node::PangraphNode;
use crate::pangraph::strand::Strand;
use itertools::Itertools;
use std::fmt::{Display, Formatter};
use std::hash::{Hash, Hasher};

#[must_use]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct SimpleNode {
  pub bid: BlockId,
  pub strand: Strand,
}

impl SimpleNode {
  pub fn new(bid: BlockId, strand: Strand) -> Self {
    Self { bid, strand }
  }

  pub fn from_full_node(n: &PangraphNode) -> Self {
    Self::new(n.block_id(), n.strand())
  }

  pub fn invert(&self) -> SimpleNode {
    SimpleNode {
      bid: self.bid,
      strand: self.strand.reverse(),
    }
  }
}

impl Display for SimpleNode {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "({}{})", self.bid, self.strand)
  }
}

#[must_use]
#[derive(Clone, Copy, Debug, Eq)]
pub struct Edge {
  pub n1: SimpleNode,
  pub n2: SimpleNode,
}

impl Edge {
  pub fn new(n1: SimpleNode, n2: SimpleNode) -> Self {
    Self { n1, n2 }
  }

  pub fn invert(&self) -> Edge {
    Edge {
      n1: self.n2.invert(),
      n2: self.n1.invert(),
    }
  }

  pub fn oriented_equal(&self, other: &Edge) -> bool {
    self.n1 == other.n1 && self.n2 == other.n2
  }
}

impl PartialEq for Edge {
  fn eq(&self, other: &Self) -> bool {
    self.oriented_equal(other) || self.oriented_equal(&other.invert())
  }
}

impl Hash for Edge {
  fn hash<H: Hasher>(&self, state: &mut H) {
    let direct_hash = {
      let mut hasher = std::collections::hash_map::DefaultHasher::new();
      (&self.n1, &self.n2).hash(&mut hasher);
      hasher.finish()
    };
    let inverted_hash = {
      let mut hasher = std::collections::hash_map::DefaultHasher::new();
      (&self.n2.invert(), &self.n1.invert()).hash(&mut hasher);
      hasher.finish()
    };
    state.write_u64(direct_hash ^ inverted_hash);
  }
}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "[{}|{}]", self.n1, self.n2)
  }
}

#[must_use]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SimplePath {
  pub nodes: Vec<SimpleNode>,
  pub circular: bool,
}

impl SimplePath {
  pub fn new(nodes: Vec<SimpleNode>, circular: bool) -> Self {
    Self { nodes, circular }
  }

  pub fn to_edges(&self) -> Vec<Edge> {
    let mut edges: Vec<Edge> = self
      .nodes
      .iter()
      .tuple_windows()
      .map(|(&n1, &n2)| Edge { n1, n2 })
      .collect();

    if self.circular {
      edges.push(Edge {
        n1: *self.nodes.last().unwrap(),
        n2: *self.nodes.first().unwrap(),
      });
    }
    edges
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use maplit::{hashmap, hashset};
  use pretty_assertions::assert_eq;

  #[test]
  fn test_simple_node() {
    let n1 = SimpleNode::new(BlockId(1), Strand::Forward);
    let n2 = SimpleNode::new(BlockId(1), Strand::Reverse);
    let n3 = SimpleNode::new(BlockId(2), Strand::Forward);
    let n4 = SimpleNode::new(BlockId(2), Strand::Reverse);

    assert_ne!(n1, n2);
    assert_eq!(n1, n2.invert());
    assert_eq!(n1.invert(), n2);
    assert_ne!(n1.invert(), n3);
    assert_ne!(n1, n3);
    assert_ne!(n1, n3.invert());
    assert_ne!(n4, n3);
    assert_eq!(n4, n3.invert());
  }

  #[test]
  fn test_edge() {
    let a = SimpleNode::new(BlockId(1), Strand::Forward);
    let b = SimpleNode::new(BlockId(2), Strand::Reverse);

    let e1 = Edge { n1: a, n2: b };

    let e2 = Edge { n1: b, n2: b };

    let e3 = Edge {
      n1: b.invert(),
      n2: a.invert(),
    };

    assert_ne!(e1, e2);
    assert_eq!(e1, e3);

    let x = hashmap! { e1 => 1, e2 => 2 };
    assert!(x.contains_key(&e3));

    assert_eq!(hashset! {e1, e2, e3}, hashset! {e1, e2});
  }
}
