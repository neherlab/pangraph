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

  use crate::graph::node::GetName;
  use crate::graph::traversal::GraphTraversalContinuation;
  use crate::io::nwk::nwk_read_str;
  use crate::io::nwk::tests::{TestEdge, TestNode};
  use eyre::Report;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::fmt::Display;
  use std::sync::Arc;

  #[test]
  fn test_postorder() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let graph = nwk_read_str::<TestNode, TestEdge>("(((A,B)AB,(C,D)CD)ABCD,(E,F)EF)root;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.traverse_backward(|node| {
      let name = node.payload.name().to_owned();
      actual.write_arc().push(name);
      GraphTraversalContinuation::Continue
    });

    assert_eq!(
      &vec!["F", "E", "D", "C", "B", "A", "EF", "CD", "AB", "ABCD", "root"],
      &*actual.read()
    );

    Ok(())
  }
}
