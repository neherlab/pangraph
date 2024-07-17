use crate::graph::edge::{EdgeFromNwk, EdgeToNwk, GraphEdge};
use crate::graph::graph::{Graph, SafeEdge, SafeNode};
use crate::graph::node::{GraphNode, GraphNodeKey, NodeFromNwk, NodeToNwk};
use crate::io::file::create_file_or_stdout;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_file_to_string;
use crate::make_error;
use crate::o;
use crate::utils::float_fmt::float_to_digits;
use bio::io::newick;
use bio_types::phylogeny::Tree;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use itertools::Itertools;
use log::warn;
use petgraph::visit::IntoNodeReferences;
use petgraph::Direction;
use smart_default::SmartDefault;
use std::io::{Cursor, Read, Write};
use std::path::Path;
use std::sync::Arc;

pub fn nwk_read_file<N, E>(filepath: impl AsRef<Path>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + NodeFromNwk,
  E: GraphEdge + EdgeFromNwk,
{
  let filepath = filepath.as_ref();
  nwk_read_reader(open_file_or_stdin(&Some(filepath))?).wrap_err_with(|| format!("When reading file '{filepath:#?}'"))
}

pub fn nwk_read_str<N, E>(s: impl AsRef<str>) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + NodeFromNwk,
  E: GraphEdge + EdgeFromNwk,
{
  let nwk_string = s.as_ref();
  nwk_read_reader(Cursor::new(nwk_string)).wrap_err_with(|| format!("When reading Newick string:\n    '{nwk_string}'"))
}

pub fn nwk_read_reader<N, E>(reader: impl Read) -> Result<Graph<N, E>, Report>
where
  N: GraphNode + NodeFromNwk,
  E: GraphEdge + EdgeFromNwk,
{
  let mut nwk_tree = newick::read(reader).wrap_err("When parsing Newick")?;

  nwk_tree.g.node_weights_mut().for_each(|weight| {
    if weight == "N/A" {
      *weight = "".to_owned();
    }
  });

  let mut graph = Graph::<N, E>::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, GraphNodeKey>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  for (nwk_idx, nwk_node) in nwk_tree.g.node_references() {
    let n_edges_incoming = nwk_tree.g.edges_directed(nwk_idx, Direction::Incoming).count();
    let n_edges_outgoing = nwk_tree.g.edges_directed(nwk_idx, Direction::Outgoing).count();

    // Discard node names which are parseable to a number. These are not names, but weights.
    // And we don't need them here. Weights are collected onto the edges later.
    let mut nwk_node: String = nwk_node.to_owned();
    if nwk_node.parse::<f64>().is_ok() {
      nwk_node = o!("");
    };

    let inserted_node_idx = match (n_edges_incoming, n_edges_outgoing) {
      (0, _) => graph.add_node(N::from_nwk(nwk_node)),
      (_, 0) => graph.add_node(N::from_nwk(nwk_node)),
      (_, _) => graph.add_node(N::from_nwk(nwk_node)),
    };

    index_map.insert(nwk_idx.index(), inserted_node_idx);
  }

  // Insert edges
  for (nwk_idx, nwk_edge) in nwk_tree.g.raw_edges().iter().enumerate() {
    let weight: f64 = nwk_edge.weight as f64;
    let source: usize = nwk_edge.source().index();
    let target: usize = nwk_edge.target().index();

    let source = index_map
      .get(&source)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {source} not found."))?;

    let target = index_map
      .get(&target)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {target} not found."))?;

    graph.add_edge(*source, *target, E::from_nwk(Some(weight)))?;
  }

  graph.build()?;

  Ok(graph)
}

#[derive(Clone, SmartDefault)]
pub struct WriteNwkOptions {
  /// Format node weights keeping this many significant digits
  pub weight_significant_digits: Option<u8>,

  /// Format node weights keeping this many decimal digits
  pub weight_decimal_digits: Option<i8>,
}

pub fn nwk_write_file<N, E>(
  filepath: &impl AsRef<Path>,
  graph: &Graph<N, E>,
  options: &WriteNwkOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
{
  let mut f = create_file_or_stdout(filepath)?;
  nwk_write_writer(&mut f, graph, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn nwk_write_str<N, E>(graph: &Graph<N, E>, options: &WriteNwkOptions) -> Result<String, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
{
  let mut buf = Vec::new();
  nwk_write_writer(&mut buf, graph, options)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nwk_write_writer<N, E>(
  writer: &mut impl Write,
  graph: &Graph<N, E>,
  options: &WriteNwkOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
{
  let roots = graph.get_roots();
  if roots.is_empty() {
    return make_error!("When converting graph to Newick format: No roots found.");
  }

  let root = {
    if roots.len() > 1 {
      unimplemented!("Multiple roots are not supported yet");
    }
    &roots[0]
  };

  let mut stack: Vec<(SafeNode<N>, Option<SafeEdge<E>>, usize)> = vec![(Arc::clone(root), None, 0)];
  while let Some((node, edge, child_visit)) = stack.pop() {
    let children = graph.children_of(&node.read()).into_iter().collect_vec();

    if child_visit < children.len() {
      stack.push((node, edge, child_visit + 1));

      if child_visit == 0 {
        write!(writer, "(")?;
      } else {
        write!(writer, ",")?;
      }

      let (child, child_edge) = &children[child_visit];
      stack.push((Arc::clone(child), Some(Arc::clone(child_edge)), 0));
    } else {
      if child_visit > 0 {
        write!(writer, ")")?;
      }

      let name = {
        let node_payload = node.read_arc().payload().read_arc();
        node_payload.to_nwk()
      };

      let weight = edge.and_then(|edge| edge.read_arc().payload().read().to_nwk());

      write!(writer, "{name}")?;

      if let Some(weight) = weight {
        write!(writer, ":{}", format_weight(weight, options))?;
      }
    }
  }

  write!(writer, ";")?;

  Ok(())
}

pub fn format_weight(weight: f64, options: &WriteNwkOptions) -> String {
  if !weight.is_finite() {
    warn!("When converting graph to Newick: Weight is invalid: '{weight}'");
  }
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),
    options.weight_decimal_digits,
  )
}

#[cfg(test)]
pub(crate) mod tests {
  use super::*;
  use crate::graph::edge::EdgeToNwk;
  use crate::graph::node::{GetName, NodeToNwk};
  use crate::io::nwk::{nwk_write_str, WriteNwkOptions};
  use derive_more::Display;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use serde::{Deserialize, Serialize};

  #[derive(Clone, Debug, Default, Display, PartialEq, Eq, Serialize, Deserialize)]
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

  #[derive(Clone, Debug, Default, Display, PartialEq, Serialize, Deserialize)]
  #[display(fmt = "{_0:?}")]
  pub struct TestEdge(pub Option<f64>);

  impl EdgeFromNwk for TestEdge {
    fn from_nwk(weight: Option<f64>) -> Self {
      Self(weight)
    }
  }

  impl EdgeToNwk for TestEdge {
    fn to_nwk(&self) -> Option<f64> {
      self.0
    }
  }

  impl GraphEdge for TestEdge {}

  #[test]
  fn test_nwk_read_write() -> Result<(), Report> {
    let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;";

    // TODO: Try to implement weights attached to root node. In a general graph, root node cannot contain
    // incoming edges. There is nowhere they could come from. So currently root node weight is not stored.
    //  Note the 0.01 at the end in this example:
    // let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    let graph = nwk_read_str::<TestNode, TestEdge>(input)?;

    let output = nwk_write_str(&graph, &WriteNwkOptions::default())?;

    assert_eq!(input, output);
    Ok(())
  }
}
