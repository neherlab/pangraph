use crate::io::fasta::{FastaReader, FastaRecord};
use eyre::{eyre, Report};
use minimap2::{Aligner, Mapping};

pub fn align_with_minimap2(reff: impl AsRef<str>, qry: impl AsRef<str>) -> Result<Vec<Vec<Mapping>>, Report> {
  let mut reader = FastaReader::from_paths(&[qry.as_ref()])?;

  let aligner = Aligner::builder()
    .asm20()
    .with_threads(8)
    .with_cigar()
    .with_index(reff.as_ref(), None)
    .map_err(|s| eyre!("Error: minimap2: when creating aligner: {s}"))?;

  let mut results = vec![];
  loop {
    let mut record = FastaRecord::default();
    reader.read(&mut record).unwrap();
    if record.is_empty() {
      break;
    }

    let alignment: Vec<Mapping> = aligner
      .map(record.seq.as_bytes(), false, false, None, None)
      .map_err(|s| {
        eyre!(
          "Error: minimap2: when aligning sequence #{} '{}': {s}",
          record.index,
          record.seq_name
        )
      })?;

    results.push(alignment);
  }

  Ok(results)
}
