use crate::commands::build::build_args::PangraphBuildArgs;
use crate::pangraph::pangraph::Pangraph;

pub fn align_self(graph: &Pangraph, args: &PangraphBuildArgs) -> Pangraph {
  // TODO: implement proper alignment
  graph.clone()
}
