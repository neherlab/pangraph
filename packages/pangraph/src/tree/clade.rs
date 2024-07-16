use crate::utils::lock::Lock;
use serde::{Deserialize, Serialize};
use std::hash::{Hash, Hasher};

#[derive(Clone, Debug, Serialize, Deserialize, Hash)]
pub struct Clade<T> {
  pub parent: Option<Lock<Clade<T>>>,
  pub left: Option<Lock<Clade<T>>>,
  pub right: Option<Lock<Clade<T>>>,
  pub data: T,
}

impl<T> Clade<T> {
  pub fn new(data: T) -> Self {
    Self {
      parent: None,
      left: None,
      right: None,
      data,
    }
  }

  pub fn from_children(data: T, left: &Lock<Clade<T>>, right: &Lock<Clade<T>>) -> Self {
    Self {
      parent: None,
      left: Some(Lock::clone(left)),
      right: Some(Lock::clone(right)),
      data,
    }
  }

  #[must_use]
  pub fn is_leaf(&self) -> bool {
    self.right.is_none() && self.left.is_none()
  }

  #[must_use]
  pub fn is_root(&self) -> bool {
    self.parent.is_none()
  }
}

pub trait WithName {
  fn name(&self) -> Option<&str>;
}

impl<T: WithName> Clade<T> {
  pub fn to_newick(&self) -> String {
    fn recurse<T: WithName>(clade: &Clade<T>) -> String {
      if clade.is_leaf() {
        String::from(clade.data.name().unwrap_or_default())
      } else {
        let mut newick = String::from("(");
        if let Some(left) = &clade.left {
          newick.push_str(&recurse(&left.read()));
        }
        newick.push(',');
        if let Some(right) = &clade.right {
          newick.push_str(&recurse(&right.read()));
        }
        newick.push(')');
        if let Some(name) = clade.data.name() {
          newick.push_str(name);
        }
        newick
      }
    }

    let newick = recurse(self);
    format!("{newick};")
  }
}

pub fn postorder<T, D, F>(clade: &Lock<Clade<D>>, f: F) -> Vec<T>
where
  F: Fn(&mut Clade<D>) -> T,
{
  fn recurse<T, D, F>(clade: &Lock<Clade<D>>, result: &mut Vec<T>, f: &F) -> ()
  where
    F: Fn(&mut Clade<D>) -> T,
  {
    if let Some(left) = &clade.read().left {
      recurse(left, result, f);
    }
    if let Some(right) = &clade.read().right {
      recurse(right, result, f);
    }
    result.push(f(&mut clade.write()));
  }

  let mut result = vec![];
  recurse(clade, &mut result, &f);
  result
}

#[cfg(test)]
mod tests {
  #![allow(clippy::many_single_char_names)]

  use super::*;
  use crate::graph::breadth_first::GraphTraversalContinuation;
  use crate::graph::create_graph_from_nwk::create_graph_from_nwk_str;
  use crate::graph::edge::{GraphEdge, Weighted};
  use crate::graph::node::{GraphNode, Named, NodeType, WithNwkComments};
  use eyre::Report;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use serde::{Deserialize, Serialize};
  use std::fmt::{Display, Formatter};
  use std::sync::Arc;

  #[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
  pub struct Node {
    pub name: String,
    pub node_type: NodeType,
  }

  impl Node {
    pub fn new(name: impl AsRef<str>, node_type: NodeType) -> Self {
      Self {
        name: name.as_ref().to_owned(),
        node_type,
      }
    }
  }

  impl GraphNode for Node {
    fn root(name: &str) -> Self {
      Self::new(name, NodeType::Root(name.to_owned()))
    }

    fn internal(name: &str) -> Self {
      Self::new(name, NodeType::Internal(name.to_owned()))
    }

    fn leaf(name: &str) -> Self {
      Self::new(name, NodeType::Leaf(name.to_owned()))
    }

    fn set_node_type(&mut self, node_type: NodeType) {
      self.node_type = node_type;
    }
  }

  impl WithNwkComments for Node {}

  impl Named for Node {
    fn name(&self) -> &str {
      &self.name
    }

    fn set_name(&mut self, name: &str) {
      self.name = name.to_owned();
    }
  }

  impl Display for Node {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
      match &self.node_type {
        NodeType::Root(weight) => write!(f, "{weight:}"),
        NodeType::Internal(weight) => write!(f, "{weight:}"),
        NodeType::Leaf(name) => write!(f, "{name}"),
      }
    }
  }

  #[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
  pub struct Edge {
    pub weight: f64,
  }

  impl GraphEdge for Edge {
    fn new(weight: f64) -> Self {
      Self { weight }
    }
  }

  impl Weighted for Edge {
    fn weight(&self) -> f64 {
      self.weight
    }
  }

  impl Display for Edge {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
      write!(f, "{:}", &self.weight)
    }
  }

  #[test]
  fn test_postorder() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let mut graph = create_graph_from_nwk_str::<Node, Edge>("(((A,B)AB,(C,D)CD)ABCD,(E,F)EF)root;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.par_iter_breadth_first_backward(|node| {
      actual.write_arc().push(node.payload.name.clone());
      GraphTraversalContinuation::Continue
    });

    assert_eq!(
      &vec!["F", "E", "D", "C", "B", "A", "EF", "CD", "AB", "ABCD", "root",],
      &*actual.read()
    );

    Ok(())
  }
}
