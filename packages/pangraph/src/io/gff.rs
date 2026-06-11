use crate::annotation::feature::{Feature, interval_from_one_based_inclusive};
use crate::io::file::open_file_or_stdin;
use crate::pangraph::strand::Strand;
use eyre::{Report, WrapErr};
use noodles::gff;
use noodles::gff::record::Strand as GffStrand;
use std::io::BufRead;
use std::path::Path;

/// Reads GFF3 annotations and normalizes each record into a format-agnostic [`Feature`].
///
/// Directives, comments and any trailing `##FASTA` section are skipped automatically.
pub struct GffReader<'a> {
  reader: gff::Reader<Box<dyn BufRead + 'a>>,
}

impl<'a> GffReader<'a> {
  /// Wrap an already-opened buffered reader.
  pub fn new(reader: Box<dyn BufRead + 'a>) -> Self {
    Self {
      reader: gff::Reader::new(reader),
    }
  }

  /// Read GFF3 from an in-memory string (handy for tests).
  pub fn from_str(contents: &'a impl AsRef<str>) -> Result<Self, Report> {
    Ok(Self::new(Box::new(contents.as_ref().as_bytes())))
  }

  /// Read GFF3 from a file path, transparently handling compression and `-` (stdin).
  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let reader = open_file_or_stdin(&Some(filepath))?;
    Ok(Self::new(reader))
  }

  /// Parse and normalize every GFF record into a [`Feature`].
  pub fn read_many(mut self) -> Result<Vec<Feature>, Report> {
    let mut features = Vec::new();
    for (i, result) in self.reader.records().enumerate() {
      let record = result.wrap_err_with(|| format!("When parsing GFF record #{}", i + 1))?;
      features.push(feature_from_gff_record(&record));
    }
    Ok(features)
  }
}

/// Convert a single noodles GFF record into our internal [`Feature`].
fn feature_from_gff_record(record: &gff::Record) -> Feature {
  let strand = match record.strand() {
    GffStrand::Forward => Some(Strand::Forward),
    GffStrand::Reverse => Some(Strand::Reverse),
    GffStrand::None | GffStrand::Unknown => None,
  };

  let attrs = record.attributes();
  let attributes: Vec<(String, String)> = attrs
    .iter()
    .map(|(tag, value)| (tag.clone(), value.to_string()))
    .collect();
  let id = attrs.get("ID").map(ToString::to_string);
  let name = attrs.get("Name").map(ToString::to_string);

  // GFF uses "." for an absent source field.
  let source = match record.source() {
    "." => None,
    s => Some(s.to_owned()),
  };

  Feature {
    seqid: record.reference_sequence_name().to_owned(),
    source,
    feature_type: record.ty().to_owned(),
    interval: interval_from_one_based_inclusive(usize::from(record.start()), usize::from(record.end())),
    strand,
    id,
    name,
    attributes,
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::utils::interval::Interval;
  use pretty_assertions::assert_eq;

  fn gff(lines: &[&str]) -> String {
    lines.join("\n")
  }

  #[test]
  fn test_gff_reader_basic_fields() -> Result<(), Report> {
    let input = gff(&[
      "##gff-version 3",
      "chr1\tprokka\tgene\t5\t30\t.\t+\t.\tID=gene0001;Name=dnaA",
    ]);
    let features = GffReader::from_str(&input)?.read_many()?;
    assert_eq!(features.len(), 1);
    let f = &features[0];
    assert_eq!(f.seqid, "chr1");
    assert_eq!(f.source, Some("prokka".to_owned()));
    assert_eq!(f.feature_type, "gene");
    assert_eq!(f.interval, Interval::new(4, 30));
    assert_eq!(f.strand, Some(Strand::Forward));
    assert_eq!(f.id, Some("gene0001".to_owned()));
    assert_eq!(f.name, Some("dnaA".to_owned()));
    assert_eq!(
      f.attributes,
      vec![
        ("ID".to_owned(), "gene0001".to_owned()),
        ("Name".to_owned(), "dnaA".to_owned())
      ]
    );
    Ok(())
  }

  #[test]
  fn test_gff_reader_strands() -> Result<(), Report> {
    let input = gff(&[
      "chr1\t.\tgene\t50\t90\t.\t-\t.\tID=g2;Name=repA",
      "chr2\t.\tmisc_feature\t10\t20\t.\t.\t.\tID=g3",
    ]);
    let features = GffReader::from_str(&input)?.read_many()?;
    assert_eq!(features.len(), 2);
    // Reverse strand, and "." source maps to None.
    assert_eq!(features[0].strand, Some(Strand::Reverse));
    assert_eq!(features[0].source, None);
    assert_eq!(features[0].name, Some("repA".to_owned()));
    // Unstranded "." maps to None, and a missing Name stays None.
    assert_eq!(features[1].seqid, "chr2");
    assert_eq!(features[1].strand, None);
    assert_eq!(features[1].id, Some("g3".to_owned()));
    assert_eq!(features[1].name, None);
    Ok(())
  }

  #[test]
  fn test_gff_reader_from_fixture_file() -> Result<(), Report> {
    let features = GffReader::from_path("../../data/example.gff")?.read_many()?;
    assert!(!features.is_empty());
    // Every feature must reference one of the two contigs in the fixture.
    assert!(features.iter().all(|f| f.seqid == "chr1" || f.seqid == "chr2"));
    Ok(())
  }
}
