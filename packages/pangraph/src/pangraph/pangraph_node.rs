use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::pangraph::strand::Strand;
use crate::utils::id::id;
use derive_more::{Display, From};
use eyre::{Context, Report};
use getset::CopyGetters;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::hash::Hash;

#[derive(
  Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize, JsonSchema,
)]
pub struct NodeId(pub usize);

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize, CopyGetters, JsonSchema)]
#[get_copy = "pub"]
pub struct PangraphNode {
  id: NodeId,
  block_id: BlockId,
  path_id: PathId,
  strand: Strand,
  position: (usize, usize),
}

impl NodeId {
  pub fn from_str(s: impl AsRef<str>) -> Result<Self, Report> {
    let s = s.as_ref();
    let id = s
      .parse::<usize>()
      .wrap_err_with(|| format!("When parsing Node ID: expected unsigned integer, but got '{s}'"))?;
    Ok(Self(id))
  }
}

impl PangraphNode {
  /// Creates a node with an explicit id.
  ///
  /// Node ids are content-derived, and [`PangraphNode::with_derived_id`] is the one place that
  /// derives them. This constructor takes the id it is given, so a caller that already has an id
  /// (because the node is being rewritten rather than created) cannot accidentally seed a second,
  /// divergent derivation scheme.
  pub fn new(id: NodeId, block_id: BlockId, path_id: PathId, strand: Strand, position: (usize, usize)) -> Self {
    Self {
      id,
      block_id,
      path_id,
      strand,
      position,
    }
  }

  /// Creates a node placing `block_id` on `path`, with an id derived from its contents.
  ///
  /// This is the only place a node id is derived. Genomes are discriminated by
  /// [`PangraphPath::seed`] rather than by their path id, so the resulting id does not depend on
  /// the order in which the input sequences were read. The path id is still what the node stores.
  pub fn with_derived_id(block_id: BlockId, path: &PangraphPath, strand: Strand, position: (usize, usize)) -> Self {
    Self {
      id: id((&block_id, &path.seed(), &strand, &position)),
      block_id,
      path_id: path.id(),
      strand,
      position,
    }
  }

  /// Moves this node onto a different path.
  ///
  /// The node id derives from the block, the genome *seed*, the strand and the position, none of
  /// which this touches, so renumbering a path leaves every node id intact. That is why this is the
  /// one field a node may have rewritten in place; everything else goes through a constructor.
  pub(crate) fn set_path_id(&mut self, path_id: PathId) {
    self.path_id = path_id;
  }

  // this is almost equivalent to checking if the node is empty
  // except for an edge case: when a circular path contains only
  // one node. In this case even if the node is not empty, the
  // start and end are the same.
  pub fn start_is_end(&self) -> bool {
    self.position.0 == self.position.1
  }
}

// #[allow(clippy::wildcard_imports)]
// mod details {
//   use super::*;
//   use eyre::eyre;
//   use serde::de::Error;
//
//   #[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
//   struct PangraphNodeWithId {
//     id: Option<NodeId>,
//     #[serde(flatten)]
//     data: PangraphNode,
//   }
//
//   impl Serialize for PangraphNode {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//       S: Serializer,
//     {
//       let dst = PangraphNodeWithId {
//         id: Some(self.id()),
//         data: self.to_owned(),
//       };
//       let s = dst.serialize(serializer)?;
//       Ok(s)
//     }
//   }
//
//   impl<'de> Deserialize<'de> for PangraphNode {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//       D: Deserializer<'de>,
//     {
//       let source = PangraphNodeWithId::deserialize(deserializer)?;
//       let id_src = source.id;
//       let id_exp = source.data.id();
//       if let Some(id_src) = id_src {
//         if id_src != id_exp {
//           return Err(Error::custom(eyre!(
//           "When deserializing node: The id property is expected to be valid hash of the node's content, but the id of the source data '{id_src}' does not match the expected id: '{id_exp}'.",
//         )));
//         }
//       } else {
//         return Err(Error::custom(eyre!(
//           "When deserializing node: The id property is expected but not found.",
//         )));
//       }
//       Ok(source.data)
//     }
//   }
// }
