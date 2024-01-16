use crate::utils::id::random_id;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
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
}

pub fn postorder<T, F>(clade: &Arc<RwLock<Clade>>, f: F) -> Vec<T>
where
  F: Fn(&Clade) -> T,
{
  let mut result = vec![];
  let mut stack = VecDeque::from([Arc::clone(clade)]);
  while let Some(current) = stack.pop_front() {
    result.push(f(&current.read_arc()));
    if let Some(right) = current.read_arc().right.as_ref() {
      stack.push_front(Arc::clone(right));
    }
    if let Some(left) = current.read_arc().left.as_ref() {
      stack.push_front(Arc::clone(left));
    }
  }
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

    let result: Vec<String> = postorder(&root, |clade| clade.name.clone().unwrap_or_default());
    assert_eq!(result, vec!["", "", "", "A", "B", "", "C", "D", "", "G", "H",]);
  }
}
