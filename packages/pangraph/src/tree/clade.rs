use crate::utils::id::random_id;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Clade {
  pub id: usize,
  pub name: Option<String>,
  pub parent: Option<Arc<RwLock<Clade>>>,
  pub left: Option<Arc<RwLock<Clade>>>,
  pub right: Option<Arc<RwLock<Clade>>>,
}

impl Clade {
  pub fn new(name: &str) -> Self {
    Self {
      id: random_id(),
      name: Some(name.to_owned()),
      parent: None,
      left: None,
      right: None,
    }
  }

  pub fn from_children(left: &Arc<RwLock<Clade>>, right: &Arc<RwLock<Clade>>) -> Self {
    Self {
      id: random_id(),
      name: None,
      parent: None,
      left: Some(Arc::clone(left)),
      right: Some(Arc::clone(right)),
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

  pub fn to_newick(&self) -> String {
    fn recurse(clade: &Clade) -> String {
      if clade.is_leaf() {
        clade.name.clone().unwrap_or_default()
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
        if let Some(name) = &clade.name {
          newick.push_str(name);
        }
        newick
      }
    }

    let newick = recurse(self);
    format!("{newick};")
  }
}

pub fn postorder<T, F>(clade: &Arc<RwLock<Clade>>, f: F) -> Vec<T>
where
  F: Fn(&Clade) -> T,
{
  fn recurse<T, F>(clade: &Arc<RwLock<Clade>>, result: &mut Vec<T>, f: &F) -> ()
  where
    F: Fn(&Clade) -> T,
  {
    if let Some(left) = &clade.read().left {
      recurse(left, result, f);
    }
    if let Some(right) = &clade.read().right {
      recurse(right, result, f);
    }
    result.push(f(&clade.read()));
  }

  let mut result = vec![];
  recurse(clade, &mut result, &f);
  result
}

#[cfg(test)]
mod tests {
  #![allow(clippy::many_single_char_names)]

  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_postorder_more_irregular_tree() {
    let a = Arc::new(RwLock::new(Clade::new("A")));
    let b = Arc::new(RwLock::new(Clade::new("B")));
    let c = Arc::new(RwLock::new(Clade::new("C")));
    let d = Arc::new(RwLock::new(Clade::new("D")));
    let e = Arc::new(RwLock::new(Clade::new("E")));
    let f = Arc::new(RwLock::new(Clade::new("F")));
    let g = Arc::new(RwLock::new(Clade::new("G")));
    let h = Arc::new(RwLock::new(Clade::new("H")));
    let i = Arc::new(RwLock::new(Clade::new("I")));
    let j = Arc::new(RwLock::new(Clade::new("J")));
    let k = Arc::new(RwLock::new(Clade::new("K")));

    let ab = Arc::new(RwLock::new(Clade::from_children(&a, &b)));
    let cd = Arc::new(RwLock::new(Clade::from_children(&c, &d)));
    let ef = Arc::new(RwLock::new(Clade::from_children(&e, &f)));
    let gh = Arc::new(RwLock::new(Clade::from_children(&g, &h)));
    let jk = Arc::new(RwLock::new(Clade::from_children(&j, &k)));

    let abcd = Arc::new(RwLock::new(Clade::from_children(&ab, &cd)));
    let root = Arc::new(RwLock::new(Clade::from_children(&abcd, &gh)));

    let nwk = root.read().to_newick();
    assert_eq!(nwk, "(((A,B),(C,D)),(G,H));");

    let result: Vec<String> = postorder(&root, |clade| clade.name.clone().unwrap_or_default());
    assert_eq!(result, vec!["A", "B", "", "C", "D", "", "", "G", "H", "", ""]);
  }
}
