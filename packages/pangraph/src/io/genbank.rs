use crate::annotation::feature::Feature;
use crate::io::file::open_file_or_stdin;
use crate::pangraph::strand::Strand;
use crate::utils::interval::Interval;
use eyre::{Report, eyre};
use gb_io::reader::SeqReader;
use gb_io::seq::{Feature as GbFeature, Location, Seq as GbSeq};
use log::warn;
use std::io::Read;
use std::path::Path;

/// Reads GenBank flat files and normalizes each feature into a format-agnostic [`Feature`].
///
/// All features of every record are emitted (including the whole-record `source` feature);
/// downstream code is free to filter by `feature_type`.
pub struct GenbankReader<'a> {
  reader: Box<dyn Read + 'a>,
}

impl<'a> GenbankReader<'a> {
  /// Wrap an already-opened reader.
  pub fn new(reader: Box<dyn Read + 'a>) -> Self {
    Self { reader }
  }

  /// Read GenBank from an in-memory string (handy for tests).
  pub fn from_str(contents: &'a impl AsRef<str>) -> Result<Self, Report> {
    Ok(Self::new(Box::new(contents.as_ref().as_bytes())))
  }

  /// Read GenBank from a file path, transparently handling compression and `-` (stdin).
  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Ok(Self::new(Box::new(open_file_or_stdin(&Some(filepath))?)))
  }

  /// Parse and normalize every feature of every record into a [`Feature`].
  pub fn read_many(self) -> Result<Vec<Feature>, Report> {
    let mut features = Vec::new();
    for result in SeqReader::new(self.reader) {
      let seq = result.map_err(|err| eyre!("When parsing GenBank record: {err}"))?;
      let seqid = genbank_seqid(&seq);
      for gb_feature in &seq.features {
        if let Some(feature) = feature_from_gb(&seqid, gb_feature) {
          features.push(feature);
        }
      }
    }
    Ok(features)
  }
}

/// The identifier used to match a GenBank record against a pangraph path: prefer the
/// `VERSION` (accession.version), then bare `ACCESSION`, then the `LOCUS` name.
fn genbank_seqid(seq: &GbSeq) -> String {
  seq
    .version
    .clone()
    .or_else(|| seq.accession.clone())
    .or_else(|| seq.name.clone())
    .unwrap_or_default()
}

/// Convert a single gb-io feature into our internal [`Feature`].
///
/// gb-io coordinates are already 0-based half-open. Compound locations (`join`/`order`)
/// are collapsed to their outer span via [`Location::find_bounds`]; a feature is on the
/// reverse strand when its location is wrapped in `complement(...)`. Features whose bounds
/// cannot be resolved are skipped with a warning.
fn feature_from_gb(seqid: &str, gb_feature: &GbFeature) -> Option<Feature> {
  let (start, end) = match gb_feature.location.find_bounds() {
    Ok(bounds) => bounds,
    Err(err) => {
      warn!(
        "Skipping GenBank feature '{}' on '{seqid}': cannot resolve location bounds: {err}",
        gb_feature.kind
      );
      return None;
    },
  };
  if start < 0 || end < start {
    warn!(
      "Skipping GenBank feature '{}' on '{seqid}': invalid bounds [{start}, {end})",
      gb_feature.kind
    );
    return None;
  }

  let strand = if matches!(gb_feature.location, Location::Complement(_)) {
    Strand::Reverse
  } else {
    Strand::Forward
  };

  let id = gb_feature.qualifier_values("locus_tag").next().map(ToOwned::to_owned);
  let name = gb_feature
    .qualifier_values("gene")
    .next()
    .or_else(|| gb_feature.qualifier_values("product").next())
    .map(ToOwned::to_owned);

  let attributes: Vec<(String, String)> = gb_feature
    .qualifiers
    .iter()
    .map(|(key, value)| (key.to_string(), value.clone().unwrap_or_default()))
    .collect();

  Some(Feature {
    seqid: seqid.to_owned(),
    source: None,
    feature_type: gb_feature.kind.to_string(),
    interval: Interval::new(start as usize, end as usize),
    strand: Some(strand),
    id,
    name,
    attributes,
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  // Build a feature-table line with the INSDC layout (key indented 5 spaces, location at
  // column 22), so the qualifier indent the parser infers is a stable 21 columns.
  fn feat(kind: &str, location: &str) -> String {
    let mut s = String::from("     ");
    s.push_str(kind);
    while s.len() < 21 {
      s.push(' ');
    }
    s.push_str(location);
    s
  }

  fn qual(q: &str) -> String {
    format!("{}{q}", " ".repeat(21))
  }

  fn example_genbank() -> String {
    let lines = vec![
      "LOCUS       chr1                     120 bp    DNA     linear   BCT 01-JAN-2020".to_owned(),
      "DEFINITION  Example record.".to_owned(),
      "ACCESSION   chr1".to_owned(),
      "VERSION     chr1.1".to_owned(),
      "FEATURES             Location/Qualifiers".to_owned(),
      feat("source", "1..120"),
      qual("/organism=\"Example\""),
      feat("gene", "5..30"),
      qual("/locus_tag=\"g0001\""),
      qual("/gene=\"dnaA\""),
      feat("CDS", "5..30"),
      qual("/locus_tag=\"g0001\""),
      qual("/gene=\"dnaA\""),
      qual("/product=\"replication initiator\""),
      feat("gene", "complement(50..90)"),
      qual("/locus_tag=\"g0002\""),
      qual("/gene=\"repA\""),
      "ORIGIN".to_owned(),
      "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt".to_owned(),
      "       61 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt".to_owned(),
      "//".to_owned(),
    ];
    let mut text = lines.join("\n");
    text.push('\n');
    text
  }

  #[test]
  fn test_genbank_seqid_prefers_version() -> Result<(), Report> {
    let input = example_genbank();
    let features = GenbankReader::from_str(&input)?.read_many()?;
    assert!(features.iter().all(|f| f.seqid == "chr1.1"));
    Ok(())
  }

  #[test]
  fn test_genbank_forward_feature_mapping() -> Result<(), Report> {
    let input = example_genbank();
    let features = GenbankReader::from_str(&input)?.read_many()?;
    let cds = features
      .iter()
      .find(|f| f.feature_type == "CDS")
      .expect("CDS feature present");
    // gb-io coordinates are already 0-based half-open: `5..30` -> [4, 30).
    assert_eq!(cds.interval, Interval::new(4, 30));
    assert_eq!(cds.strand, Some(Strand::Forward));
    assert_eq!(cds.id, Some("g0001".to_owned()));
    assert_eq!(cds.name, Some("dnaA".to_owned()));
    assert!(
      cds
        .attributes
        .contains(&("product".to_owned(), "replication initiator".to_owned()))
    );
    Ok(())
  }

  #[test]
  fn test_genbank_complement_is_reverse() -> Result<(), Report> {
    let input = example_genbank();
    let features = GenbankReader::from_str(&input)?.read_many()?;
    let rep = features
      .iter()
      .find(|f| f.id.as_deref() == Some("g0002"))
      .expect("repA feature present");
    assert_eq!(rep.interval, Interval::new(49, 90));
    assert_eq!(rep.strand, Some(Strand::Reverse));
    assert_eq!(rep.name, Some("repA".to_owned()));
    Ok(())
  }

  #[test]
  fn test_genbank_reader_from_fixture_file() -> Result<(), Report> {
    let features = GenbankReader::from_path("../../data/example.gbk")?.read_many()?;
    assert!(!features.is_empty());
    assert!(features.iter().all(|f| f.seqid == "chr1.1"));
    // The fixture's repA gene is on the complement strand.
    let rep = features
      .iter()
      .find(|f| f.id.as_deref() == Some("g0002"))
      .expect("repA feature present");
    assert_eq!(rep.strand, Some(Strand::Reverse));
    Ok(())
  }
}
