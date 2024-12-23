use crate::io::file::create_file_or_stdout;
use crate::io::write::WriteAdapterIoToFmt;
use crate::o;
use crate::pangraph::pangraph::Pangraph;
use eyre::{Context, Report};
use gfa::gfa::Orientation::{Backward, Forward};
use gfa::gfa::{Link, Path as GFAPath, Segment, GFA};
use gfa::optfields::OptionalFields;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

pub fn gfa_write_file(filepath: impl AsRef<Path>, g: &Pangraph) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  gfa_write(create_file_or_stdout(filepath)?, g).wrap_err_with(|| format!("When writing gfa file: {filepath:#?}"))
}

pub fn gfa_write_str(g: &Pangraph) -> Result<String, Report> {
  let mut buf = vec![];
  gfa_write(&mut buf, g)?;
  Ok(String::from_utf8(buf)?)
}

pub fn gfa_write<W: Write>(writer: W, g: &Pangraph) -> Result<(), Report> {
  let gfa = convert_pangraph_to_gfa(g)?;
  gfa::writer::write_gfa(&gfa, &mut WriteAdapterIoToFmt(writer));
  Ok(())
}

fn convert_pangraph_to_gfa(pangraph: &Pangraph) -> Result<GFA<Vec<u8>, OptionalFields>, Report> {
  let mut gfa = GFA::<Vec<u8>, OptionalFields>::new();

  gfa.segments = pangraph
    .blocks
    .par_iter()
    .map(|(_, block)| Segment {
      name: block.id().to_string().as_bytes().to_vec(),
      sequence: block.consensus().as_bytes().to_vec(),
      optional: OptionalFields::new(),
    })
    .collect();

  gfa.links = pangraph
    .nodes
    .par_iter()
    .map(|(_, node)| {
      // FIXME: bogus data
      Link::<Vec<u8>, OptionalFields> {
        from_segment: o!("A").as_bytes().to_vec(),
        from_orient: Backward,
        to_segment: o!("B").as_bytes().to_vec(),
        to_orient: Forward,
        overlap: b"1M".to_vec(),
        optional: OptionalFields::new(),
      }
    })
    .collect();

  gfa.paths = pangraph
    .paths
    .par_iter()
    .map(|(_, path)| {
      let segment_names = path
        .nodes
        .iter()
        .map(|node_id| pangraph.nodes[node_id].block_id().to_string())
        .collect::<Vec<_>>();

      let path_name = path
        .name
        .clone()
        .unwrap_or_else(|| path.id().to_string())
        .as_bytes()
        .to_vec();

      let segment_names = segment_names.join(",").as_bytes().to_vec();
      let overlaps = vec![None; segment_names.len() - 1];
      let optional = OptionalFields::new();

      GFAPath::new(path_name, segment_names, overlaps, optional)
    })
    .collect();

  Ok(gfa)
}

#[cfg(test)]
mod tests {
  use super::*;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;

  #[test]
  fn test_gfa_empty() {
    let actual = gfa_write_str(&Pangraph::default()).unwrap();
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

    let actual = gfa_write_str(&g).unwrap();

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

    assert_eq!(actual, expected);
  }
}
