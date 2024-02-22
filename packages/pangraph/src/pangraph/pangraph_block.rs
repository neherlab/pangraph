use crate::align::alignment::Alignment;
use crate::align::map_variations::map_variations;
use crate::pangraph::edits::Edits;
use crate::utils::id::random_id;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct PangraphBlock {
  pub id: usize,
  pub consensus: String,
  pub alignments: BTreeMap<usize, Edits>,
}

impl PangraphBlock {
  pub fn from_consensus(consensus: impl AsRef<str>) -> Self {
    Self {
      id: random_id(),
      consensus: consensus.as_ref().to_owned(),
      alignments: BTreeMap::new(),
    }
  }

  pub fn new(consensus: impl AsRef<str>, alignments: BTreeMap<usize, Edits>) -> Self {
    Self {
      id: random_id(),
      consensus: consensus.as_ref().to_owned(),
      alignments,
    }
  }

  /// Append a sequence to a block, and assign the given node_id
  pub fn append_sequence(&mut self, new_seq: impl AsRef<str>, new_node_id: usize) -> Result<(), Report> {
    // Calculate the set of edits
    let edits = map_variations(&self.consensus, new_seq)?;
    // Append the edits to the block alignment dictionary
    self.alignments.insert(new_node_id, edits);
    Ok(())
  }
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::pangraph::edits::{Del, Ins, Sub};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[test]
  fn test_block_append_sequence_simple_case() {
    let consensus = o!("ACTTTGCGTATTTACTATA");

    let E_0 = Edits {
      dels: vec![Del::new(2, 3)],
      subs: vec![],
      inss: vec![],
    };
    let E_1 = Edits {
      inss: vec![Ins::new(16, "CCC")],
      subs: vec![],
      dels: vec![],
    };
    let E_2 = Edits {
      subs: vec![Sub::new(11, 'C')],
      dels: vec![],
      inss: vec![],
    };

    let alignments = BTreeMap::from([(0, E_0), (1, E_1), (2, E_2)]);

    let mut block = PangraphBlock {
      id: 0,
      consensus,
      alignments,
    };

    let new_seq = o!("ACTTTGCGGATTTACTATA");
    block.append_sequence(new_seq, 0).unwrap();

    let E_new_expected = Edits {
      subs: vec![Sub::new(8, 'G')],
      dels: vec![],
      inss: vec![],
    };

    assert_eq!(block.alignments[&0], E_new_expected);
  }
}
