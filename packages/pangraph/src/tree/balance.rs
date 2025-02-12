use crate::tree::clade::Clade;
use crate::utils::lock::Lock;

// Rotate a binary tree to a balanced configuration.
// Preserves the topological ordering of the leaves.
pub fn balance<T: Default>(clade: &Lock<Clade<T>>) -> Lock<Clade<T>> {
  let tips = leaves(clade);
  bisect(&tips)
}

fn bisect<T: Default>(tips: &[Lock<Clade<T>>]) -> Lock<Clade<T>> {
  if tips.len() <= 1 {
    return Lock::clone(&tips[0]);
  }
  let m = tips.len() / 2;
  let l = bisect(&tips[..m]);
  let r = bisect(&tips[m..]);
  Lock::new(Clade::from_children(T::default(), &l, &r))
}

fn leaves<T>(root: &Lock<Clade<T>>) -> Vec<Lock<Clade<T>>> {
  fn recurse<T>(node: &Lock<Clade<T>>, result: &mut Vec<Lock<Clade<T>>>) {
    if node.read().is_leaf() {
      result.push(node.clone());
    }

    if let Some(left) = &node.read().left {
      recurse(left, result);
    }

    if let Some(right) = &node.read().right {
      recurse(right, result);
    }
  }

  let mut result = vec![];
  recurse(root, &mut result);
  result
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::tree::clade::WithName;
  use eyre::Report;
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

  #[allow(clippy::many_single_char_names)]
  #[rstest]
  fn test_balance() -> Result<(), Report> {
    //               __ A
    //      ________|
    //     |        |__ H
    //     |
    //     |          __ B
    //     |       __|
    //   __|      |  |__ E
    //  |  |   ___|
    //  |  |  |   |__ D
    //  |  |__|
    //  |     |    __ C
    //  |     |___|
    //  |         |__ G
    //  |
    //  |____________ F

    let root = {
      let a = Lock::new(Clade::new(N::new("A")));
      let b = Lock::new(Clade::new(N::new("B")));
      let c = Lock::new(Clade::new(N::new("C")));
      let d = Lock::new(Clade::new(N::new("D")));
      let e = Lock::new(Clade::new(N::new("E")));
      let f = Lock::new(Clade::new(N::new("F")));
      let g = Lock::new(Clade::new(N::new("G")));
      let h = Lock::new(Clade::new(N::new("H")));
      let i = Lock::new(Clade::new(N::new("I")));
      let j = Lock::new(Clade::new(N::new("J")));
      let k = Lock::new(Clade::new(N::new("K")));

      let ah = Lock::new(Clade::from_children(N::new(""), &a, &h));
      let be = Lock::new(Clade::from_children(N::new(""), &b, &e));
      let bed = Lock::new(Clade::from_children(N::new(""), &be, &d));
      let cg = Lock::new(Clade::from_children(N::new(""), &c, &g));
      let bedcg = Lock::new(Clade::from_children(N::new(""), &bed, &cg));
      let ahbedcg = Lock::new(Clade::from_children(N::new(""), &ah, &bedcg));
      Lock::new(Clade::from_children(N::new(""), &ahbedcg, &f))
    };

    let actual = balance(&root);
    let actual = actual.read().to_newick();
    assert_eq!(actual, "(((A,H),(B,E)),((D,C),(G,F)));");

    Ok(())
  }
}
