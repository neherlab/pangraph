use crate::graph::node::GraphNodeKey;
use derive_more::Display;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::hash::Hash;
use std::sync::Arc;

pub trait GraphEdge: Clone + Debug + Sync + Send {}

pub trait EdgeFromNwk {
  fn from_nwk(weight: Option<f64>) -> Self;
}

pub trait EdgeToNwk {
  fn to_nwk(&self) -> Option<f64>;
}

pub trait GetWeight {
  fn weight(&self) -> Option<f64>;
}

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct GraphEdgeKey(pub usize);

impl GraphEdgeKey {
  #[inline]
  pub const fn as_usize(self) -> usize {
    self.0
  }
}

/// Edge representing a connection between two nodes. Relevant data can be
/// stored in the edge atomically. Edge's target and source node's are
/// weak references and can't outlive the nodes they represent.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Edge<E: GraphEdge> {
  key: GraphEdgeKey,
  source: GraphNodeKey,
  target: GraphNodeKey,
  data: Arc<RwLock<E>>,
}

impl<E: GraphEdge> Edge<E> {
  /// Creates a new edge.
  pub fn new(key: GraphEdgeKey, source: GraphNodeKey, target: GraphNodeKey, data: E) -> Edge<E> {
    Edge {
      key,
      source,
      target,
      data: Arc::new(RwLock::new(data)),
    }
  }

  #[inline]
  pub const fn key(&self) -> GraphEdgeKey {
    self.key
  }

  #[inline]
  pub const fn source(&self) -> GraphNodeKey {
    self.source
  }

  #[inline]
  pub const fn target(&self) -> GraphNodeKey {
    self.target
  }

  #[inline]
  pub fn payload(&self) -> Arc<RwLock<E>> {
    Arc::clone(&self.data)
  }
}
