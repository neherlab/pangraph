use crate::utils::lock::Lock;
use serde::{Deserialize, Serialize};
use std::hash::Hash;

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
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[derive(Default)]
  struct N(pub String);

  impl N {
    pub fn new(name: impl Into<String>) -> Self {
      Self(name.into())
    }
  }

  impl WithName for N {
    fn name(&self) -> Option<&str> {
      Some(&self.0)
    }
  }

  #[rstest]
  fn test_postorder_more_irregular_tree() {
    let a = Lock::new(Clade::new(N::new("A")));
    let b = Lock::new(Clade::new(N::new("B")));
    let c = Lock::new(Clade::new(N::new("C")));
    let d = Lock::new(Clade::new(N::new("D")));
    let _e = Lock::new(Clade::new(N::new("E")));
    let _f = Lock::new(Clade::new(N::new("F")));
    let g = Lock::new(Clade::new(N::new("G")));
    let h = Lock::new(Clade::new(N::new("H")));
    let _i = Lock::new(Clade::new(N::new("I")));
    let _j = Lock::new(Clade::new(N::new("J")));
    let _k = Lock::new(Clade::new(N::new("K")));

    let ab = Lock::new(Clade::from_children(N::new(""), &a, &b));
    let cd = Lock::new(Clade::from_children(N::new(""), &c, &d));
    let gh = Lock::new(Clade::from_children(N::new(""), &g, &h));

    let abcd = Lock::new(Clade::from_children(N::new(""), &ab, &cd));
    let root = Lock::new(Clade::from_children(N::new(""), &abcd, &gh));

    let nwk = root.read().to_newick();
    assert_eq!(nwk, "(((A,B),(C,D)),(G,H));");

    let result: Vec<String> = postorder(&root, |clade| clade.data.0.clone());
    assert_eq!(result, vec!["A", "B", "", "C", "D", "", "", "G", "H", "", ""]);
  }
}
