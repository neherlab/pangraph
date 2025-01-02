use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_path::PathId;
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
  pub fn new(
    node_id: Option<NodeId>,
    block_id: BlockId,
    path_id: PathId,
    strand: Strand,
    position: (usize, usize),
  ) -> Self {
    let id = node_id.unwrap_or_else(|| id((&block_id, &path_id, &strand, &position)));
    Self {
      id,
      block_id,
      path_id,
      strand,
      position,
    }
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
