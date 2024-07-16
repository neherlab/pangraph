use crate::io::fasta::FastaRecord;
use crate::make_internal_error;
use crate::pangraph::pangraph::Pangraph;
use eyre::Report;

pub fn verify_result(graph: &Pangraph, fastas: &[FastaRecord]) -> Result<(), Report> {
  if graph.paths.len() != fastas.len() {
    return make_internal_error!(
      "Resulting graph contains {} paths, but there were {} input sequences",
      graph.paths.len(),
      fastas.len()
    );
  }

  // TODO: implement the actual verification logic

  Ok(())
}
