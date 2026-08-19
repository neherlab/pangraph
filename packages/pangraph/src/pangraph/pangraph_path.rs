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

  #[getset(get = "pub")]
  pub desc: Option<String>,
}

impl PangraphPath {
  /// Creates a path with an explicit id.
  ///
  /// Path ids are sequential rather than content-derived: they double as the ordering index of the
  /// genomes, so `build` assigns them from the input order and `renumber_paths` shifts them during
  /// a merge. The id is therefore always the caller's to supply, and there is no fallback that
  /// could seed a second, divergent numbering.
  pub fn new(
    id: PathId,
    nodes: impl Into<Vec<NodeId>>,
    tot_len: usize,
    circular: bool,
    name: Option<String>,
    desc: Option<String>,
  ) -> Self {
    Self {
      id,
      nodes: nodes.into(),
      tot_len,
      circular,
      name,
      desc,
    }
  }

  /// Order-independent identity of the genome on this path, used as the discriminator when deriving
  /// node ids.
  ///
  /// Taken from the genome name rather than from the path id, so that node ids do not depend on the
  /// order in which the input sequences were read, and do not change when `renumber_paths` shifts
  /// path ids during a merge. Genome names are unique within a graph, and `merge` rejects graphs
  /// that share one, so the seed identifies a genome as well as the path id does.
  ///
  /// Recomputed on demand rather than stored: a stored field would have to be kept out of the JSON
  /// and then recomputed on load anyway, and silently defaults to the same value for every path if
  /// that is forgotten. Falls back to the path id for unnamed paths, which pangraph never writes.
  pub fn seed(&self) -> usize {
    self.name.as_deref().map_or(self.id.0, genome_seed)
  }
}

/// Derives the identifier seed of a genome from its name.
///
/// The single definition of the seeding rule, so that [`PangraphPath::seed`] and
/// `Pangraph::singleton`, which needs the seed before it has a path to ask, cannot drift apart.
pub fn genome_seed(name: &str) -> usize {
  id(name)
}
