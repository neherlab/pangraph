use crate::commands::build::build_args::PangraphBuildArgs;
use crate::pangraph::pangraph::Pangraph;

pub fn align_pair(left: &Pangraph, right: &Pangraph, args: &PangraphBuildArgs) -> Pangraph {
  // TODO: implement proper alignment
  left.clone()
}
