use crate::pangraph::pangraph_node::NodeId;
use crate::utils::id::id;
use derive_more::{Display, From};
use getset::{CopyGetters, Getters};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(
  Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize, JsonSchema,
)]
pub struct PathId(pub usize);

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, Getters, CopyGetters, JsonSchema)]
pub struct PangraphPath {
  #[getset(get_copy = "pub")]
  pub id: PathId,

  #[getset(get = "pub")]
  pub nodes: Vec<NodeId>,

  #[getset(get_copy = "pub")]
  pub tot_len: usize,

  #[getset(get_copy = "pub")]
  pub circular: bool,

  #[getset(get = "pub")]
  pub name: Option<String>,
}

impl PangraphPath {
  pub fn new(
    path_id: Option<PathId>,
    nodes: impl Into<Vec<NodeId>>,
    tot_len: usize,
    circular: bool,
    name: Option<String>,
  ) -> Self {
    let nodes = nodes.into();
    let id = path_id.unwrap_or_else(|| id((&nodes, &tot_len, &circular)));
    Self {
      id,
      nodes,
      tot_len,
      circular,
      name,
    }
  }
}
