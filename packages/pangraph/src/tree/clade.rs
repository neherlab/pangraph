use crate::utils::id::random_id;
use crate::utils::lock::Lock;
use serde::{Deserialize, Serialize};
use std::hash::{Hash, Hasher};

#[derive(Clone, Debug, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct Clade {
  pub id: usize,
  pub name: Option<String>,
  pub parent: Option<Lock<Clade>>,
  pub left: Option<Lock<Clade>>,
  pub right: Option<Lock<Clade>>,
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

  pub fn from_children(left: &Lock<Clade>, right: &Lock<Clade>) -> Self {
    Self {
      id: random_id(),
      name: None,
      parent: None,
      left: Some(Lock::clone(left)),
      right: Some(Lock::clone(right)),
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

pub fn postorder<T, F>(clade: &Lock<Clade>, f: F) -> Vec<T>
where
  F: Fn(&Clade) -> T,
{
  fn recurse<T, F>(clade: &Lock<Clade>, result: &mut Vec<T>, f: &F) -> ()
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
    let a = Lock::new(Clade::new("A"));
    let b = Lock::new(Clade::new("B"));
    let c = Lock::new(Clade::new("C"));
    let d = Lock::new(Clade::new("D"));
    let e = Lock::new(Clade::new("E"));
    let f = Lock::new(Clade::new("F"));
    let g = Lock::new(Clade::new("G"));
    let h = Lock::new(Clade::new("H"));
    let i = Lock::new(Clade::new("I"));
    let j = Lock::new(Clade::new("J"));
    let k = Lock::new(Clade::new("K"));

    let ab = Lock::new(Clade::from_children(&a, &b));
    let cd = Lock::new(Clade::from_children(&c, &d));
    let gh = Lock::new(Clade::from_children(&g, &h));

    let abcd = Lock::new(Clade::from_children(&ab, &cd));
    let root = Lock::new(Clade::from_children(&abcd, &gh));

    let nwk = root.read().to_newick();
    assert_eq!(nwk, "(((A,B),(C,D)),(G,H));");

    let result: Vec<String> = postorder(&root, |clade| clade.name.clone().unwrap_or_default());
    assert_eq!(result, vec!["A", "B", "", "C", "D", "", "", "G", "H", "", ""]);
  }
}
