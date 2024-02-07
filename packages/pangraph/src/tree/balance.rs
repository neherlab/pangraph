use crate::tree::clade::Clade;
use crate::utils::lock::Lock;
use itertools::Itertools;

// Rotate a binary tree to a balanced configuration.
// Preserves the topological ordering of the leaves.
pub fn balance(clade: &Lock<Clade>) -> Lock<Clade> {
  let tips = leaves(clade);
  bisect(&tips)
}

fn bisect(tips: &[Lock<Clade>]) -> Lock<Clade> {
  if tips.len() <= 1 {
    return Lock::clone(&tips[0]);
  }
  let m = tips.len() / 2;
  let l = bisect(&tips[..m]);
  let r = bisect(&tips[m..]);
  Lock::new(Clade::from_children(&l, &r))
}

fn leaves(root: &Lock<Clade>) -> Vec<Lock<Clade>> {
  fn recurse(node: &Lock<Clade>, result: &mut Vec<Lock<Clade>>) {
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
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

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

      let ah = Lock::new(Clade::from_children(&a, &h));
      let be = Lock::new(Clade::from_children(&b, &e));
      let bed = Lock::new(Clade::from_children(&be, &d));
      let cg = Lock::new(Clade::from_children(&c, &g));
      let bedcg = Lock::new(Clade::from_children(&bed, &cg));
      let ahbedcg = Lock::new(Clade::from_children(&ah, &bedcg));
      Lock::new(Clade::from_children(&ahbedcg, &f))
    };

    let actual = balance(&root);
    let actual = actual.read().to_newick();
    assert_eq!(actual, "(((A,H),(B,E)),((D,C),(G,F)));");

    Ok(())
  }
}
