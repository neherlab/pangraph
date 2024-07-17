use crate::graph::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph::node::{GraphNode, GraphNodeKey, Node};
use crate::graph::traversal::{
  directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward, GraphTraversalContinuation,
};
use crate::{make_error, make_internal_error, make_internal_report};
use eyre::Report;
use itertools::Itertools;
use parking_lot::lock_api::{ArcRwLockReadGuard, ArcRwLockWriteGuard};
use parking_lot::{RawRwLock, RwLock};
use serde::{Deserialize, Serialize};
use std::fmt::{Debug, Display};
use std::io::Write;
use std::sync::Arc;

pub type SafeNode<N> = Arc<RwLock<Node<N>>>;
pub type SafeNodeRef<N> = ArcRwLockReadGuard<RawRwLock, Node<N>>;
pub type SafeNodeRefMut<N> = ArcRwLockWriteGuard<RawRwLock, Node<N>>;

pub type SafeEdge<E> = Arc<RwLock<Edge<E>>>;
pub type SafeEdgeRef<E> = ArcRwLockReadGuard<RawRwLock, Edge<E>>;
pub type SafeEdgeRefMut<E> = ArcRwLockWriteGuard<RawRwLock, Edge<E>>;

pub type SafeNodePayloadRef<N> = ArcRwLockReadGuard<RawRwLock, N>;
pub type SafeNodePayloadRefMut<N> = ArcRwLockWriteGuard<RawRwLock, N>;

pub type NodeEdgePair<N, E> = (Arc<RwLock<Node<N>>>, Arc<RwLock<Edge<E>>>);
pub type NodeEdgePayloadPair<N, E> = (Arc<RwLock<N>>, Arc<RwLock<E>>);

/// Represents graph node during forward traversal
#[derive(Debug)]
pub struct GraphNodeForward<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: SafeNodePayloadRefMut<N>,
  pub parents: Vec<NodeEdgePayloadPair<N, E>>,
}

impl<N, E> GraphNodeForward<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub fn new(graph: &Graph<N, E>, node: &Node<N>) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let payload = node.payload().write_arc();

    let parents = graph
      .parents_of(node)
      .iter()
      .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
      .collect_vec();

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      parents,
    }
  }
}

/// Represents graph node during backwards traversal
#[derive(Debug)]
pub struct GraphNodeBackward<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: SafeNodePayloadRefMut<N>,
  pub children: Vec<NodeEdgePayloadPair<N, E>>,
}

impl<N, E> GraphNodeBackward<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub fn new(graph: &Graph<N, E>, node: &Node<N>) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let payload = node.payload().write_arc();

    let children = graph
      .children_of(node)
      .iter()
      .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
      .collect_vec();

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
    }
  }
}

/// Represents graph node during safe traversal
#[derive(Debug)]
pub struct GraphNodeSafe<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: Arc<RwLock<N>>,
  pub children: Vec<NodeEdgePair<N, E>>,
  pub parents: Vec<NodeEdgePair<N, E>>,
}

impl<N, E> GraphNodeSafe<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub fn from_node(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Self {
    let node = node.read();
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();
    let payload = node.payload();
    let parents = graph.parents_of(&node);
    let children = graph.children_of(&node);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      parents,
    }
  }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Graph<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  nodes: Vec<Arc<RwLock<Node<N>>>>,
  edges: Vec<Arc<RwLock<Edge<E>>>>,
  roots: Vec<GraphNodeKey>,
  leaves: Vec<GraphNodeKey>,
}

impl<N, E> Graph<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub const fn new() -> Self {
    Self {
      nodes: Vec::new(),
      edges: Vec::new(),
      roots: vec![],
      leaves: vec![],
    }
  }

  /// Retrieve parent nodes of a given node and the corresponding edges.
  ///
  /// **Returns**: list of pairs `(parent, edge)`, where `parent` is the parent node,
  /// and `edge` is the inbound edge connecting the parent node with the given node.
  pub fn parents_of(&self, node: &Node<N>) -> Vec<NodeEdgePair<N, E>> {
    // Parents are the source nodes of inbound edges
    node
      .inbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .filter_map(|edge| {
        let parent_key = edge.read().source();
        self.get_node(parent_key).map(|parent| (parent, edge))
      })
      .collect_vec()
  }

  /// Retrieve keys of parent nodes of a given node.
  pub fn parent_keys_of(&self, node: &Node<N>) -> Vec<GraphNodeKey> {
    // Parents are the source nodes of inbound edges
    node
      .inbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .map(|edge| edge.read().source())
      .collect_vec()
  }

  pub fn exactly_one_parent_of(&self, node: &Node<N>) -> Result<Arc<RwLock<Node<N>>>, Report> {
    self.one_parent_of(node)?.ok_or_else(|| {
      make_internal_report!(
        "No parents found for node {} (context: is_root={} is_leaf={})",
        node.key(),
        node.is_root(),
        node.is_leaf()
      )
    })
  }

  pub fn one_parent_of(&self, node: &Node<N>) -> Result<Option<Arc<RwLock<Node<N>>>>, Report> {
    let parents = self
      .parents_of(node)
      .into_iter()
      .map(|(parent, _)| parent)
      .collect_vec();

    if parents.is_empty() {
      return Ok(None);
    }

    if parents.len() > 1 {
      return make_internal_error!(
        "Only trees with exactly one parent per node are currently supported, but node '{}' has {} parents",
        node.key(),
        self.roots.len()
      );
    }

    Ok(Some(Arc::clone(&parents[0])))
  }

  /// Retrieve child nodes of a given node and the corresponding edges.
  ///
  /// **Returns**: list of pairs `(child, edge)`, where `child` is the child node,
  /// and `edge` is the outbound edge connecting the given node with the child node.
  pub fn children_of(&self, node: &Node<N>) -> Vec<NodeEdgePair<N, E>> {
    // Children are the target nodes of outbound edges
    node
      .outbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .filter_map(|edge| {
        let child_key = edge.read().target();
        self.get_node(child_key).map(|child| (child, edge))
      })
      .collect_vec()
  }

  /// Retrieve keys of child nodes of a given node.
  pub fn child_keys_of(&self, node: &Node<N>) -> Vec<GraphNodeKey> {
    // Children are the target nodes of outbound edges
    node
      .outbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .map(|edge| edge.read().target())
      .collect_vec()
  }

  /// Retrieve keys of child nodes of a given node by key
  pub fn child_keys_of_by_key(&self, node_key: GraphNodeKey) -> Vec<GraphNodeKey> {
    let node = self.get_node(node_key).unwrap().read_arc();
    self.child_keys_of(&node)
  }

  pub fn get_node(&self, index: GraphNodeKey) -> Option<Arc<RwLock<Node<N>>>> {
    self.nodes.get(index.0).map(Arc::clone)
  }

  pub fn get_edge(&self, index: GraphEdgeKey) -> Option<Arc<RwLock<Edge<E>>>> {
    self.edges.get(index.0).map(Arc::clone)
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn for_each<T, F>(&self, f: &mut dyn FnMut(GraphNodeSafe<N, E>)) {
    self
      .nodes
      .iter()
      .for_each(|node| f(GraphNodeSafe::from_node(self, node)));
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E>) -> T,
  {
    self
      .nodes
      .iter()
      .map(|node| f(GraphNodeSafe::from_node(self, node)))
      .collect_vec()
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn filter_map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E>) -> Option<T>,
  {
    self
      .nodes
      .iter()
      .filter_map(|node| f(GraphNodeSafe::from_node(self, node)))
      .collect_vec()
  }

  #[inline]
  pub fn num_nodes(&self) -> usize {
    self.nodes.len()
  }

  #[inline]
  pub fn num_roots(&self) -> usize {
    self.roots.len()
  }

  #[inline]
  pub fn num_leaves(&self) -> usize {
    self.leaves.len()
  }

  #[inline]
  pub fn get_nodes(&self) -> &[Arc<RwLock<Node<N>>>] {
    &self.nodes
  }

  #[inline]
  pub fn get_node_payloads(&self) -> impl Iterator<Item = Arc<RwLock<N>>> + '_ {
    self.nodes.iter().map(|node| node.read().payload())
  }

  pub fn get_exactly_one_root(&self) -> Result<Arc<RwLock<Node<N>>>, Report> {
    let roots = self.get_roots();
    if roots.len() != 1 {
      make_internal_error!(
        "Only trees with exactly one root are currently supported, but found '{}'",
        self.roots.len()
      )
    } else {
      Ok(Arc::clone(&roots[0]))
    }
  }

  #[inline]
  pub fn get_roots(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec()
  }

  #[inline]
  pub fn get_root_payloads(&self) -> impl Iterator<Item = Arc<RwLock<N>>> + '_ {
    self
      .roots
      .iter()
      .filter_map(|idx| self.get_node(*idx))
      .map(|node| node.read().payload())
  }

  #[inline]
  pub fn get_leaves(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self.leaves.iter().filter_map(|idx| self.get_node(*idx)).collect_vec()
  }

  #[inline]
  pub fn get_edges(&self) -> &[Arc<RwLock<Edge<E>>>] {
    &self.edges
  }

  pub fn add_node(&mut self, node_payload: N) -> GraphNodeKey {
    let node_key = GraphNodeKey(self.nodes.len());
    let node = Arc::new(RwLock::new(Node::new(node_key, node_payload)));
    self.nodes.push(node);
    node_key
  }

  /// Add a new edge to the graph.
  pub fn add_edge(
    &mut self,
    source_key: GraphNodeKey,
    target_key: GraphNodeKey,
    edge_payload: E,
  ) -> Result<(), Report> {
    if source_key == target_key {
      return make_error!(
        "When adding a graph edge {source_key}->{target_key}: Attempted to connect node {source_key} to itself."
      );
    }

    let source_lock = self
      .get_node(source_key)
      .ok_or_else(|| format!("When adding a graph edge {source_key}->{target_key}: Node {source_key} not found."))
      .unwrap();

    let target_lock = self
      .get_node(target_key)
      .ok_or_else(|| format!("When adding a graph edge {source_key}->{target_key}: Node {target_key} not found."))
      .unwrap();

    let edge_key = GraphEdgeKey(self.edges.len());
    let new_edge = Arc::new(RwLock::new(Edge::new(edge_key, source_key, target_key, edge_payload)));

    {
      let (source, target) = (source_lock.read(), target_lock.read());

      let already_connected = source
        .outbound()
        .iter()
        .any(|edge| self.get_edge(*edge).unwrap().read().target() == target.key());

      if already_connected {
        return make_error!("When adding a graph edge {source_key}->{target_key}: Nodes {source_key} and {target_key} are already connected.");
      }

      self.edges.push(Arc::clone(&new_edge));
    }

    {
      let (mut source, mut target) = (source_lock.write(), target_lock.write());
      source.outbound_mut().push(edge_key);
      target.inbound_mut().push(edge_key);
    }

    Ok(())
  }

  pub fn build(&mut self) -> Result<(), Report> {
    self.roots = self
      .nodes
      .iter()
      .filter_map(|node| {
        let node = node.read();
        node.is_root().then(|| node.key())
      })
      .collect_vec();

    self.leaves = self
      .nodes
      .iter()
      .filter_map(|node| {
        let node = node.read();
        node.is_leaf().then(|| node.key())
      })
      .collect_vec();

    Ok(())
  }

  pub fn traverse_forward<F>(&self, explorer: F)
  where
    F: Fn(GraphNodeForward<N, E>) -> GraphTraversalContinuation + Sync + Send,
  {
    let roots = self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec();

    directed_breadth_first_traversal_forward::<N, E, _>(self, roots.as_slice(), |node| {
      explorer(GraphNodeForward::new(self, node))
    });

    self.reset_nodes();
  }

  pub fn traverse_backward<F>(&self, explorer: F)
  where
    F: Fn(GraphNodeBackward<N, E>) -> GraphTraversalContinuation + Sync + Send,
  {
    let leaves = self
      .leaves
      .iter()
      .filter_map(|idx| self.get_node(*idx))
      .rev()
      .collect_vec();

    directed_breadth_first_traversal_backward::<N, E, _>(self, leaves.as_slice(), |node| {
      explorer(GraphNodeBackward::new(self, node))
    });

    self.reset_nodes();
  }

  /// Returns graph into initial state after traversal
  pub fn reset_nodes(&self) {
    // Mark all nodes as not visited. As a part of traversal all nodes, one by one, marked as visited,
    // to ensure correctness of traversal. Here we reset the "is visited" markers this,
    // to allow for traversals again.
    self.nodes.iter().for_each(|node| node.write().mark_as_not_visited());
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::edge::{EdgeFromNwk, EdgeToNwk};
  use crate::graph::node::{GetName, NodeFromNwk, NodeToNwk};
  use crate::io::nwk::nwk_read_str;
  use derive_more::Display;

  #[derive(Clone, Debug, Display, PartialEq, Eq, Serialize, Deserialize)]
  pub struct TestNode(pub String);

  impl GraphNode for TestNode {}

  impl NodeToNwk for TestNode {
    fn to_nwk(&self) -> String {
      self.0.clone()
    }
  }

  impl NodeFromNwk for TestNode {
    fn from_nwk(name: String) -> Self {
      Self(name)
    }
  }

  impl GetName for TestNode {
    fn name(&self) -> &str {
      &self.0
    }
  }

  #[derive(Clone, Debug, Display, Eq, PartialEq, Serialize, Deserialize)]
  pub struct TestEdge;

  impl EdgeFromNwk for TestEdge {
    fn from_nwk(weight: Option<f64>) -> Self {
      Self
    }
  }

  impl EdgeToNwk for TestEdge {
    fn to_nwk(&self) -> Option<f64> {
      None
    }
  }

  impl GraphEdge for TestEdge {}

  #[test]
  fn test_traversal_parallel_breadth_first_forward() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let graph = nwk_read_str::<TestNode, TestEdge>("(((A,B)AB,(C,D)CD)ABCD,(E,F)EF)root;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.traverse_forward(|node| {
      actual.write_arc().push(node.payload.name().to_owned());
      GraphTraversalContinuation::Continue
    });

    assert_eq!(
      &vec!["root", "ABCD", "EF", "AB", "CD", "E", "F", "A", "B", "C", "D"],
      &*actual.read()
    );

    Ok(())
  }

  #[test]
  fn test_traversal_parallel_breadth_first_backward() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let graph = nwk_read_str::<TestNode, TestEdge>("(((A,B)AB,(C,D)CD)ABCD,(E,F)EF)root;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.traverse_backward(|node| {
      actual.write_arc().push(node.payload.name().to_owned());
      GraphTraversalContinuation::Continue
    });

    assert_eq!(
      &vec!["F", "E", "D", "C", "B", "A", "EF", "CD", "AB", "ABCD", "root"],
      &*actual.read()
    );

    Ok(())
  }
}
