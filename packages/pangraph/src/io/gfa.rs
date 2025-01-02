use crate::io::file::create_file_or_stdout;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::strand::Strand;
use clap::Parser;
use derive_more::Constructor;
use eyre::{Context, Report};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::Path;

#[derive(Parser, Debug, Default, Clone)]
pub struct GfaWriteParams {
  /// Blocks below this length cutoff will not be exported
  #[clap(long)]
  pub minimum_length: Option<usize>,

  /// Blocks above this length cutoff will not be exported
  #[clap(long)]
  pub maximum_length: Option<usize>,

  /// Blocks below this depth cutoff will not be exported
  #[clap(long)]
  pub minimum_depth: Option<usize>,

  /// Blocks above this depth cutoff will not be exported
  #[clap(long)]
  pub maximum_depth: Option<usize>,

  /// Include block sequences in the GFA file
  #[clap(long)]
  pub include_sequences: bool,

  /// Exclude blocks that are duplicated in any path
  #[clap(long)]
  pub no_duplicated: bool,
}

impl GfaWriteParams {
  pub fn get_minimum_length(&self) -> usize {
    self.minimum_length.unwrap_or(0)
  }

  pub fn get_maximum_length(&self) -> usize {
    self.maximum_length.unwrap_or(usize::MAX)
  }

  pub fn get_minimum_depth(&self) -> usize {
    self.minimum_depth.unwrap_or(0)
  }

  pub fn get_maximum_depth(&self) -> usize {
    self.maximum_depth.unwrap_or(usize::MAX)
  }

  pub fn include_sequences(&self) -> bool {
    self.include_sequences
  }

  pub fn no_duplicated(&self) -> bool {
    self.no_duplicated
  }
}

pub fn gfa_write_file(filepath: impl AsRef<Path>, g: &Pangraph, params: &GfaWriteParams) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  gfa_write(create_file_or_stdout(filepath)?, g, params)
    .wrap_err_with(|| format!("When writing gfa file: {filepath:#?}"))
}

pub fn gfa_write_str(g: &Pangraph, params: &GfaWriteParams) -> Result<String, Report> {
  let mut buf = vec![];
  gfa_write(&mut buf, g, params)?;
  Ok(String::from_utf8(buf)?)
}

pub fn gfa_write<W: Write>(mut writer: W, g: &Pangraph, params: &GfaWriteParams) -> Result<(), Report> {
  let gfa = convert_pangraph_to_gfa(g, params)?;

  writeln!(writer, "H\tVN:Z:1.0")?;

  if !gfa.segments.is_empty() {
    writeln!(writer, "# blocks")?;
  }

  for segment in gfa.segments.values() {
    let segment_seq = if params.include_sequences() {
      &segment.sequence
    } else {
      "*"
    };
    let name = segment.name;
    let read_count = segment.depth * segment.length;
    let length = segment.length;
    let duplicated_tag = if segment.duplicated { "\tDP:Z:duplicated" } else { "" };
    writeln!(
      writer,
      "S\t{name}\t{segment_seq}\tRC:i:{read_count}\tLN:i:{length}{duplicated_tag}",
    )?;
  }

  if !gfa.links.edge_ct.is_empty() {
    writeln!(writer, "# edges")?;
  }

  for (edge, read_count) in &gfa.links.edge_ct {
    let bid1 = edge.n1.bid;
    let strand1 = edge.n1.strand;
    let bid2 = edge.n2.bid;
    let strand2 = edge.n2.strand;
    writeln!(writer, "L\t{bid1}\t{strand1}\t{bid2}\t{strand2}\t*\tRC:i:{read_count}",)?;
  }

  if !gfa.paths.is_empty() {
    writeln!(writer, "# paths")?;
  }
  for path in &gfa.paths {
    let circular_tag = if path.circular { "\tTP:Z:circular" } else { "" };
    let segments = path
      .segments
      .iter()
      .map(|node| format!("{}{}", node.bid, node.strand))
      .join(",");
    let path_name = &path.path_name;
    writeln!(writer, "P\t{path_name}\t{segments}\t*{circular_tag}",)?;
  }
  Ok(())
}

#[derive(Debug, Clone, Default)]
pub struct Gfa {
  segments: BTreeMap<BlockId, GfaSegment>,
  links: GfaLinks,
  paths: Vec<GfaPath>,
}

#[derive(Debug, Clone)]
pub struct GfaSegment {
  name: BlockId,
  sequence: String,
  depth: usize,
  length: usize,
  duplicated: bool,
}

#[derive(Debug, Clone, Hash, Eq, PartialEq, Ord, PartialOrd, Constructor)]
pub struct SimpleNode {
  bid: BlockId,
  strand: Strand,
}

#[derive(Debug, Clone, Hash, Eq, PartialEq, Ord, PartialOrd, Constructor)]
pub struct Edge {
  n1: SimpleNode,
  n2: SimpleNode,
}

#[derive(Debug, Clone, Default)]
pub struct GfaLinks {
  edge_ct: BTreeMap<Edge, usize>,
}

impl GfaLinks {
  pub fn add_edge(&mut self, bid1: BlockId, strand1: Strand, bid2: BlockId, strand2: Strand) {
    let n1 = SimpleNode::new(bid1, strand1);
    let n2 = SimpleNode::new(bid2, strand2);
    let edge = Edge::new(n1, n2);
    *self.edge_ct.entry(edge).or_insert(0) += 1;
  }
}

#[derive(Debug, Clone)]
pub struct GfaPath {
  path_name: String,
  segments: Vec<SimpleNode>,
  circular: bool,
}

fn convert_pangraph_to_gfa(graph: &Pangraph, params: &GfaWriteParams) -> Result<Gfa, Report> {
  let segments = graph_to_gfa_segment(graph);

  let paths = graph_to_gfa_paths(graph)
    .iter()
    .map(|path| gfa_filter_path(path, &segments, params))
    .filter(|path| !path.segments.is_empty())
    .collect_vec();

  let segment_ids = paths
    .iter()
    .flat_map(|path| path.segments.iter().map(|segment| segment.bid))
    .collect::<BTreeSet<_>>();

  let segments = segments.into_iter().filter(|(k, _)| segment_ids.contains(k)).collect();

  let gfa_links = gfa_paths_to_links(&paths);

  Ok(Gfa {
    segments,
    links: gfa_links,
    paths,
  })
}

fn gfa_paths_to_links(paths: &[GfaPath]) -> GfaLinks {
  let mut links = GfaLinks::default();
  for path in paths {
    for window in path.segments.windows(2) {
      assert!(window.len() > 1);
      let (seg1, seg2) = (&window[0], &window[1]);
      links.add_edge(seg1.bid, seg1.strand, seg2.bid, seg2.strand);
    }
    if path.circular {
      let seg1 = path.segments.last().unwrap();
      let seg2 = &path.segments[0];
      links.add_edge(seg1.bid, seg1.strand, seg2.bid, seg2.strand);
    }
  }
  links
}

fn graph_to_gfa_paths(graph: &Pangraph) -> Vec<GfaPath> {
  let mut paths = Vec::new();
  for path in graph.paths.values() {
    let segments = path
      .nodes
      .iter()
      .map(|node_id| {
        let node = &graph.nodes[node_id];
        SimpleNode::new(node.block_id(), node.strand())
      })
      .collect();

    paths.push(GfaPath {
      path_name: path.name.as_ref().unwrap().clone(),
      segments,
      circular: path.circular,
    });
  }
  paths
}

fn gfa_filter_path(path: &GfaPath, segments: &BTreeMap<BlockId, GfaSegment>, params: &GfaWriteParams) -> GfaPath {
  let new_segments: Vec<_> = path
    .segments
    .iter()
    .filter(|node| {
      let segment = &segments[&node.bid];
      let length = segment.length;
      let depth = segment.depth;

      length >= params.get_minimum_length()
        && length <= params.get_maximum_length()
        && depth >= params.get_minimum_depth()
        && depth <= params.get_maximum_depth()
        && (!params.no_duplicated() || !segment.duplicated)
    })
    .cloned()
    .collect();

  GfaPath {
    path_name: path.path_name.clone(),
    segments: new_segments,
    circular: path.circular,
  }
}

fn graph_to_gfa_segment(graph: &Pangraph) -> BTreeMap<BlockId, GfaSegment> {
  graph
    .blocks
    .iter()
    .map(|(block_id, block)| {
      let segment = GfaSegment {
        name: *block_id,
        sequence: block.consensus().to_owned(),
        depth: block.depth(),
        length: block.consensus_len(),
        duplicated: block.is_duplicated(graph),
      };
      (*block_id, segment)
    })
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;

  #[test]
  fn test_gfa_empty() {
    let actual = gfa_write_str(&Pangraph::default(), &GfaWriteParams::default()).unwrap();
    let expected = indoc! {r#"
    H	VN:Z:1.0
    "#};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_gfa_general_case() {
    let g = Pangraph::from_str(indoc! {
    // language=json
    r#"
    {
      "paths": {
        "0": {
          "id": 0,
          "nodes": [
            14515840915932838377
          ],
          "tot_len": 1737,
          "circular": false,
          "name": "Path A"
        },
        "1": {
          "id": 1,
          "nodes": [
            15291847754458130853
          ],
          "tot_len": 1737,
          "circular": false,
          "name": "Path B"
        },
        "2": {
          "id": 2,
          "nodes": [
            15109482180931348145
          ],
          "tot_len": 1737,
          "circular": false,
          "name": "Path C"
        }
      },
      "blocks": {
        "12778560093473594666": {
          "id": 12778560093473594666,
          "consensus": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
          "alignments": {
            "14515840915932838377": {
              "subs": [],
              "dels": [],
              "inss": []
            },
            "15291847754458130853": {
              "subs": [],
              "dels": [],
              "inss": []
            },
            "15109482180931348145": {
              "subs": [],
              "dels": [],
              "inss": []
            }
          }
        }
      },
      "nodes": {
        "14515840915932838377": {
          "id": 14515840915932838377,
          "block_id": 12778560093473594666,
          "path_id": 0,
          "strand": "+",
          "position": [
            0,
            0
          ]
        },
        "15291847754458130853": {
          "id": 15291847754458130853,
          "block_id": 12778560093473594666,
          "path_id": 1,
          "strand": "+",
          "position": [
            0,
            0
          ]
        },
        "15109482180931348145": {
          "id": 15109482180931348145,
          "block_id": 12778560093473594666,
          "path_id": 2,
          "strand": "+",
          "position": [
            0,
            0
          ]
        }
      }
    }
    "#})
    .unwrap();

    let params = GfaWriteParams {
      include_sequences: true,
      ..Default::default()
    };

    let actual = gfa_write_str(&g, &params).unwrap();

    let expected = indoc! {r#"
    H	VN:Z:1.0
    S	12778560093473594666	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    L	A	-	B	+	1M
    L	A	-	B	+	1M
    L	A	-	B	+	1M
    P	Path A	12778560093473594666	*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*
    P	Path B	12778560093473594666	*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*
    P	Path C	12778560093473594666	*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*
    "#};

    assert_eq!(expected, actual);
  }
}
