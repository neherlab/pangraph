use crate::commands::reconstruct::reconstruct_args::PangraphReconstructArgs;
use crate::io::fasta::{FastaReader, FastaRecord, FastaWriter};
use crate::io::json::json_read_file;
use crate::make_error;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::reconstruct::reconstruct;
use eyre::Report;
use log::info;

pub fn reconstruct_run(args: &PangraphReconstructArgs) -> Result<(), Report> {
  let PangraphReconstructArgs {
    input_graph,
    output_fasta,
    verify,
  } = &args;

  let graph: Pangraph = json_read_file(input_graph)?;
  let mut results = reconstruct(&graph);

  if let Some(verify) = verify {
    info!("Verifying sequences reconstructed from pangenome graph");
    let mut reader = FastaReader::from_path(verify)?;
    results.try_for_each(|actual| -> Result<(), Report> {
      let actual = actual?;
      let mut expected = FastaRecord::new();
      reader.read(&mut expected)?;
      compare_sequences(&expected, &actual)?;
      Ok(())
    })?;
  } else {
    let mut writer = FastaWriter::from_path(output_fasta)?;
    results.try_for_each(|fasta| {
      let fasta = fasta?;
      writer.write(fasta.seq_name, &fasta.desc, &fasta.seq)
    })?;
  }

  Ok(())
}

pub fn compare_sequences(left: &FastaRecord, right: &FastaRecord) -> Result<bool, Report> {
  if left != right {
    return make_error!(
      "Sequence mismatch detected: expected length {} but got {}",
      left.seq.len(),
      right.seq.len()
    );
  }
  Ok(true)
}
