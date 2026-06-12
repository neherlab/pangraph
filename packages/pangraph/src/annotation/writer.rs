use crate::annotation::lift::LiftedAnnotation;
use crate::io::file::create_file_or_stdout;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_node::NodeId;
use crate::pangraph::strand::Strand;
use csv::{Writer as CsvWriter, WriterBuilder};
use eyre::Report;
use serde::Serialize;
use std::io::Write;
use std::path::Path;

/// A pluggable sink for lifted annotations.
///
/// Concrete implementations choose the output format (CSV is the default; a nested JSON writer
/// arrives with the block-level phase) without the lift logic knowing anything about it — the
/// lift only ever produces [`LiftedAnnotation`] objects.
pub trait AnnotationWriter {
  /// Serialize a batch of node-level lifted annotations.
  fn write_node_annotations(&mut self, annotations: &[LiftedAnnotation]) -> Result<(), Report>;
}

/// One CSV row per [`LiftedAnnotation`].
///
/// Borrows from the source annotation to avoid cloning. `attributes` is rendered as a single
/// JSON-string column (a JSON array of `[key, value]` pairs, order- and duplicate-preserving);
/// `Option` fields render as empty cells; strong-typed ids render as their plain numbers;
/// `strand_on_consensus` renders as `+`/`-`/empty.
#[derive(Serialize)]
struct LiftedAnnotationCsvRow<'a> {
  feature_id: &'a str,
  parent_feature_id: Option<&'a str>,
  segment_idx: usize,
  n_segments: usize,
  genome: &'a str,
  block_id: BlockId,
  node_id: NodeId,
  strand_on_consensus: Option<Strand>,
  node_start: usize,
  node_end: usize,
  cons_start: usize,
  cons_end: usize,
  start_is_terminus: bool,
  end_is_terminus: bool,
  start_in_insertion: bool,
  end_in_insertion: bool,
  frac_covered: String,
  #[serde(rename = "type")]
  feature_type: &'a str,
  name: Option<&'a str>,
  attributes: String,
}

impl<'a> LiftedAnnotationCsvRow<'a> {
  fn from_lifted(a: &'a LiftedAnnotation) -> Result<Self, Report> {
    Ok(Self {
      feature_id: &a.feature_id,
      parent_feature_id: a.parent_feature_id.as_deref(),
      segment_idx: a.segment_idx,
      n_segments: a.n_segments,
      genome: &a.genome,
      block_id: a.block_id,
      node_id: a.node_id,
      strand_on_consensus: a.strand_on_consensus,
      node_start: a.node_start,
      node_end: a.node_end,
      cons_start: a.cons_start,
      cons_end: a.cons_end,
      start_is_terminus: a.start_is_terminus,
      end_is_terminus: a.end_is_terminus,
      start_in_insertion: a.start_in_insertion,
      end_in_insertion: a.end_in_insertion,
      frac_covered: format!("{:.4}", a.frac_covered),
      feature_type: &a.feature_type,
      name: a.name.as_deref(),
      attributes: serde_json::to_string(&a.attributes)?,
    })
  }
}

/// The default [`AnnotationWriter`]: long-format CSV, one row per lifted segment.
///
/// Backed by [`create_file_or_stdout`], so `-` writes to stdout and the output is transparently
/// compressed when the path carries a `.gz`/`.bz2`/`.xz`/`.zst` extension. The header row is
/// written on the first record.
pub struct CsvAnnotationWriter {
  writer: CsvWriter<Box<dyn Write + Send>>,
}

impl CsvAnnotationWriter {
  /// Create a CSV annotation writer at `filepath` (`-` = stdout). `delimiter` is typically `b','`.
  pub fn new(filepath: impl AsRef<Path>, delimiter: u8) -> Result<Self, Report> {
    let file = create_file_or_stdout(filepath)?;
    let writer = WriterBuilder::new()
      .delimiter(delimiter)
      .has_headers(true)
      .from_writer(file);
    Ok(Self { writer })
  }
}

impl AnnotationWriter for CsvAnnotationWriter {
  fn write_node_annotations(&mut self, annotations: &[LiftedAnnotation]) -> Result<(), Report> {
    for ann in annotations {
      let row = LiftedAnnotationCsvRow::from_lifted(ann)?;
      self.writer.serialize(&row)?;
    }
    self.writer.flush()?;
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::io::csv::parse_csv;
  use serde::Deserialize;
  use std::fs::read_to_string;
  use tempfile::tempdir;

  /// Minimal owned mirror of the CSV row for round-tripping the written file back in.
  /// Some fields are only used to drive deserialization, not asserted on.
  #[allow(dead_code)]
  #[derive(Debug, Deserialize, PartialEq)]
  struct Row {
    feature_id: String,
    parent_feature_id: Option<String>,
    segment_idx: usize,
    n_segments: usize,
    genome: String,
    block_id: usize,
    node_id: usize,
    strand_on_consensus: Option<String>,
    node_start: usize,
    node_end: usize,
    cons_start: usize,
    cons_end: usize,
    start_is_terminus: bool,
    end_is_terminus: bool,
    start_in_insertion: bool,
    end_in_insertion: bool,
    frac_covered: String,
    #[serde(rename = "type")]
    feature_type: String,
    name: Option<String>,
    attributes: String,
  }

  fn sample(feature_id: &str, segment_idx: usize) -> LiftedAnnotation {
    LiftedAnnotation {
      feature_id: feature_id.to_owned(),
      parent_feature_id: Some("g1".to_owned()),
      segment_idx,
      n_segments: 2,
      genome: "genomeA".to_owned(),
      block_id: BlockId(7),
      node_id: NodeId(42),
      strand_on_consensus: Some(Strand::Reverse),
      node_start: 3,
      node_end: 8,
      cons_start: 3,
      cons_end: 9,
      start_is_terminus: true,
      end_is_terminus: false,
      start_in_insertion: false,
      end_in_insertion: true,
      frac_covered: 0.5,
      feature_type: "CDS".to_owned(),
      name: Some("geneA".to_owned()),
      attributes: vec![
        ("ID".to_owned(), "g1".to_owned()),
        ("Name".to_owned(), "geneA".to_owned()),
      ],
    }
  }

  #[test]
  fn test_csv_writer_round_trips_rows() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("annotations.csv");

    let anns = vec![sample("g1.seg0", 0), sample("g1.seg1", 1)];
    {
      let mut writer = CsvAnnotationWriter::new(&path, b',').unwrap();
      writer.write_node_annotations(&anns).unwrap();
    }

    let contents = read_to_string(&path).unwrap();
    // Header is written, and the attributes column is a JSON string.
    assert!(contents.starts_with("feature_id,parent_feature_id,segment_idx"));
    assert!(contents.contains(r#"[[""ID"",""g1""],[""Name"",""geneA""]]"#));

    let rows: Vec<Row> = parse_csv(&contents).unwrap();
    assert_eq!(rows.len(), 2);
    let r = &rows[0];
    assert_eq!(r.feature_id, "g1.seg0");
    assert_eq!(r.parent_feature_id.as_deref(), Some("g1"));
    assert_eq!((r.block_id, r.node_id), (7, 42));
    assert_eq!(r.strand_on_consensus.as_deref(), Some("-"));
    assert_eq!((r.cons_start, r.cons_end), (3, 9));
    assert_eq!((r.start_is_terminus, r.end_is_terminus), (true, false));
    assert_eq!((r.start_in_insertion, r.end_in_insertion), (false, true));
    assert_eq!(r.frac_covered, "0.5000");
    assert_eq!(r.feature_type, "CDS");
    assert_eq!(r.attributes, r#"[["ID","g1"],["Name","geneA"]]"#);
  }

  #[test]
  fn test_csv_writer_empty_input_writes_no_header() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("empty.csv");
    {
      let mut writer = CsvAnnotationWriter::new(&path, b',').unwrap();
      writer.write_node_annotations(&[]).unwrap();
    }
    assert_eq!(read_to_string(&path).unwrap(), "");
  }
}
