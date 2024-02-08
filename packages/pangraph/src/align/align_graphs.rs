use crate::align::align_pair::align_pair;
use crate::align::align_self::align_self;
use crate::commands::build::build_args::PangraphBuildArgs;
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn align_graphs(left: &Pangraph, right: &Pangraph, args: &PangraphBuildArgs) -> Result<Pangraph, Report> {
  let graph = align_pair(left, right, args);
  let graph = align_self(&graph, args);
  Ok(graph)
}
