use crate::commands::reconstruct::reconstruct_args::PangraphReconstructArgs;
use crate::io::fasta::{FastaReader, FastaRecord, FastaWriter};
use crate::io::json::json_read_file;
use crate::make_error;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::reconstruct::{path_ids_by_name, reconstruct, reconstruct_genome, verify_genome};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;
use std::collections::BTreeSet;
use std::path::Path;

pub fn reconstruct_run(args: &PangraphReconstructArgs) -> Result<(), Report> {
  let PangraphReconstructArgs {
    input_graph,
    output_fasta,
    verify,
  } = &args;

  let graph: Pangraph = json_read_file(input_graph)?;

  if let Some(verify) = verify {
    info!("Verifying sequences reconstructed from pangenome graph");
    let n_verified = verify_against_fasta(&graph, verify)
      .wrap_err_with(|| format!("When verifying reconstructed sequences against '{}'", verify.display()))?;
    info!(
      "Graph reconstructs all {n_verified} genomes of '{}' exactly",
      verify.display()
    );
  } else {
    let mut writer = FastaWriter::from_path(output_fasta)?;
    reconstruct(&graph).try_for_each(|fasta| {
      let fasta = fasta?;
      writer.write(fasta.seq_name, &fasta.desc, &fasta.seq)
    })?;
  }

  Ok(())
}

/// Checks that the graph reconstructs exactly the genomes of the given FASTA file, and returns how
/// many were verified.
///
/// Genomes are matched by name, so the order of the FASTA records is irrelevant: a merged graph
/// renumbers its path ids and reconstructs its genomes in an order that matches no input file.
/// The file is streamed and genomes are reconstructed one at a time, so only a single genome is
/// held in memory at once.
fn verify_against_fasta(graph: &Pangraph, verify: &Path) -> Result<usize, Report> {
  let path_ids = path_ids_by_name(graph)?;
  let mut remaining: BTreeSet<&str> = path_ids.keys().copied().collect();

  let mut reader = FastaReader::from_path(verify)?;
  let mut record = FastaRecord::new();
  loop {
    // `FastaReader::read` clears the record before filling it, so it is not reset here.
    reader.read(&mut record)?;
    if record.is_empty() {
      break;
    }

    let Some(path_id) = path_ids.get(record.seq_name.as_str()) else {
      return make_error!(
        "Verification file contains genome '{}', which the graph does not contain",
        record.seq_name
      );
    };

    if !remaining.remove(record.seq_name.as_str()) {
      return make_error!(
        "Verification file contains genome '{}' more than once. Genome names must be unique, because they identify genomes.",
        record.seq_name
      );
    }

    let actual = reconstruct_genome(graph, *path_id)
      .wrap_err_with(|| format!("When reconstructing genome '{}'", record.seq_name))?;
    verify_genome(&record.seq_name, &record.seq, &actual)?;
  }

  if !remaining.is_empty() {
    return make_error!(
      "Verification file is missing {} genome(s) that the graph contains: [{}]",
      remaining.len(),
      remaining.iter().join(", ")
    );
  }

  Ok(path_ids.len())
}
