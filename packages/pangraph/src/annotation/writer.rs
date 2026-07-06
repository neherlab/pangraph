use crate::annotation::compact::BlockAnnotation;
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

  /// Serialize a batch of block-level (compacted) annotations.
  fn write_block_annotations(&mut self, annotations: &[BlockAnnotation]) -> Result<(), Report>;
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
  feature_strand: Option<Strand>,
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
      feature_strand: a.feature_strand,
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

/// One CSV row per [`BlockAnnotation`] — i.e. one segment of a compacted crossing.
///
/// Mirrors [`LiftedAnnotationCsvRow`]: borrows from the source, renders `consensus_attributes` as a
/// single JSON-string column (a JSON array of `[key, value]` pairs), `strand_on_consensus` as
/// `+`/`-`/empty, ids as plain numbers, and `Option`s as empty cells. Rows of one crossing share a
/// `cluster_id` and are ordered 5'→3' by `segment_idx`.
#[derive(Serialize)]
struct BlockAnnotationCsvRow<'a> {
  #[serde(rename = "type")]
  feature_type: &'a str,
  cluster_id: usize,
  segment_idx: usize,
  n_segments: usize,
  block_id: BlockId,
  cons_start: usize,
  cons_end: usize,
  strand_on_consensus: Option<Strand>,
  consensus_name: Option<&'a str>,
  consensus_attributes: String,
  n_support: usize,
  n_total: usize,
  n_support_segment: usize,
  n_total_segment: usize,
}

impl<'a> BlockAnnotationCsvRow<'a> {
  fn from_block(a: &'a BlockAnnotation) -> Result<Self, Report> {
    Ok(Self {
      feature_type: &a.feature_type,
      cluster_id: a.cluster_id,
      segment_idx: a.segment_idx,
      n_segments: a.n_segments,
      block_id: a.block_id,
      cons_start: a.cons_start,
      cons_end: a.cons_end,
      strand_on_consensus: a.strand_on_consensus,
      consensus_name: a.consensus_name.as_deref(),
      consensus_attributes: serde_json::to_string(&a.consensus_attributes)?,
      n_support: a.n_support,
      n_total: a.n_total,
      n_support_segment: a.n_support_segment,
      n_total_segment: a.n_total_segment,
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

  fn write_block_annotations(&mut self, annotations: &[BlockAnnotation]) -> Result<(), Report> {
    for ann in annotations {
      let row = BlockAnnotationCsvRow::from_block(ann)?;
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
  use crate::pangraph::pangraph_path::PathId;
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
    feature_strand: Option<String>,
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
      path_id: PathId(2),
      block_id: BlockId(7),
      node_id: NodeId(42),
      strand_on_consensus: Some(Strand::Reverse),
      feature_strand: Some(Strand::Forward),
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
    assert_eq!(r.feature_strand.as_deref(), Some("+"));
    assert_eq!((r.cons_start, r.cons_end), (3, 9));
    assert_eq!((r.start_is_terminus, r.end_is_terminus), (true, false));
    assert_eq!((r.start_in_insertion, r.end_in_insertion), (false, true));
    assert_eq!(r.frac_covered, "0.5000");
    assert_eq!(r.feature_type, "CDS");
    assert_eq!(r.attributes, r#"[["ID","g1"],["Name","geneA"]]"#);
  }

  /// Minimal owned mirror of the block-level CSV row.
  #[allow(dead_code)]
  #[derive(Debug, Deserialize, PartialEq)]
  struct BlockRow {
    #[serde(rename = "type")]
    feature_type: String,
    cluster_id: usize,
    segment_idx: usize,
    n_segments: usize,
    block_id: usize,
    cons_start: usize,
    cons_end: usize,
    strand_on_consensus: Option<String>,
    consensus_name: Option<String>,
    consensus_attributes: String,
    n_support: usize,
    n_total: usize,
    n_support_segment: usize,
    n_total_segment: usize,
  }

  fn sample_block() -> BlockAnnotation {
    BlockAnnotation {
      feature_type: "CDS".to_owned(),
      cluster_id: 4,
      segment_idx: 1,
      n_segments: 2,
      block_id: BlockId(7),
      cons_start: 3,
      cons_end: 120,
      strand_on_consensus: Some(Strand::Reverse),
      consensus_name: Some("geneA".to_owned()),
      consensus_attributes: vec![("product".to_owned(), "widget".to_owned())],
      n_support: 18,
      n_total: 21,
      n_support_segment: 19,
      n_total_segment: 25,
    }
  }

  #[test]
  fn test_csv_writer_round_trips_block_rows() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("block_annotations.csv");

    let anns = vec![sample_block()];
    {
      let mut writer = CsvAnnotationWriter::new(&path, b',').unwrap();
      writer.write_block_annotations(&anns).unwrap();
    }

    let contents = read_to_string(&path).unwrap();
    assert!(contents.starts_with("type,cluster_id,segment_idx"));
    assert!(contents.contains(r#"[[""product"",""widget""]]"#));

    let rows: Vec<BlockRow> = parse_csv(&contents).unwrap();
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(r.feature_type, "CDS");
    assert_eq!((r.cluster_id, r.segment_idx, r.n_segments), (4, 1, 2));
    assert_eq!(r.block_id, 7);
    assert_eq!((r.cons_start, r.cons_end), (3, 120));
    assert_eq!(r.strand_on_consensus.as_deref(), Some("-"));
    assert_eq!(r.consensus_name.as_deref(), Some("geneA"));
    assert_eq!((r.n_support, r.n_total), (18, 21));
    assert_eq!((r.n_support_segment, r.n_total_segment), (19, 25));
    assert_eq!(r.consensus_attributes, r#"[["product","widget"]]"#);
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
