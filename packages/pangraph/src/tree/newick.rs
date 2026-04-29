use crate::io::file::open_file_or_stdin;
use crate::make_error;
use crate::pangraph::pangraph::Pangraph;
use crate::tree::clade::{Clade, WithName};
use crate::utils::lock::Lock;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;

impl<T: WithName> Clade<T> {
  pub fn to_newick(&self) -> String {
    fn recurse<T: WithName>(clade: &Clade<T>) -> String {
      if clade.is_leaf() {
        String::from(clade.data.name().unwrap_or_default())
      } else {
        let mut newick = String::from("(");
        if let Some(left) = &clade.left {
          newick.push_str(&recurse(&left.read()));
        }
        newick.push(',');
        if let Some(right) = &clade.right {
          newick.push_str(&recurse(&right.read()));
        }
        newick.push(')');
        if let Some(name) = clade.data.name() {
          newick.push_str(name);
        }
        newick
      }
    }

    let newick = recurse(self);
    format!("{newick};")
  }
}

/// Parse Newick text into a tree of leaf names. Internal nodes carry `None`; leaves carry
/// `Some(name)`. Branch lengths and internal labels are accepted but ignored. Only strictly
/// bifurcating trees are accepted.
pub fn parse_newick(s: &str) -> Result<Lock<Clade<Option<String>>>, Report> {
  let mut parser = Parser::new(s);
  parser.skip_ws();
  if parser.eof() {
    return make_error!("Newick input is empty");
  }
  let root = parser.parse_subtree()?;
  parser.skip_ws();
  if parser.peek() == Some(';') {
    parser.advance();
  }
  parser.skip_ws();
  if !parser.eof() {
    return make_error!(
      "Newick: unexpected trailing content at position {}: '{}'",
      parser.pos,
      parser.remaining()
    );
  }
  Ok(root)
}

/// Read a Newick file (transparent decompression), parse it, and attach singleton graphs
/// to its leaves by matching FASTA `seq_name` against leaf labels.
///
/// Validates that the tree is strictly bifurcating and that every input genome appears
/// exactly once as a leaf.
pub fn build_tree_from_newick<P: AsRef<Path>>(
  path: P,
  graphs: Vec<Pangraph>,
) -> Result<Lock<Clade<Option<Pangraph>>>, Report> {
  let path = path.as_ref();
  let path_opt: Option<&Path> = Some(path);
  let mut reader =
    open_file_or_stdin(&path_opt).wrap_err_with(|| format!("When opening guide tree file '{}'", path.display()))?;
  let mut text = String::new();
  reader
    .read_to_string(&mut text)
    .wrap_err_with(|| format!("When reading guide tree file '{}'", path.display()))?;
  let name_tree = parse_newick(&text).wrap_err_with(|| format!("When parsing Newick from '{}'", path.display()))?;

  let mut by_name: BTreeMap<String, Pangraph> = BTreeMap::new();
  for graph in graphs {
    let name = singleton_seq_name(&graph)?;
    if by_name.insert(name.clone(), graph).is_some() {
      return make_error!("Duplicate FASTA sequence name '{name}'");
    }
  }

  let tree = attach_graphs(&name_tree, &mut by_name)?;

  if !by_name.is_empty() {
    let leftover = by_name.keys().cloned().collect::<Vec<_>>().join(", ");
    return make_error!("FASTA records [{leftover}] are not present in the guide tree");
  }

  Ok(tree)
}

fn singleton_seq_name(graph: &Pangraph) -> Result<String, Report> {
  graph
    .paths
    .values()
    .next()
    .and_then(|p| p.name.clone())
    .ok_or_else(|| eyre::eyre!("Encountered a singleton graph without a path name"))
}

fn attach_graphs(
  src: &Lock<Clade<Option<String>>>,
  by_name: &mut BTreeMap<String, Pangraph>,
) -> Result<Lock<Clade<Option<Pangraph>>>, Report> {
  let (left_opt, right_opt, name_opt) = {
    let g = src.read();
    (g.left.clone(), g.right.clone(), g.data.clone())
  };
  match (left_opt, right_opt) {
    (None, None) => {
      let leaf_name = name_opt.ok_or_else(|| eyre::eyre!("Newick leaf without a name"))?;
      let graph = by_name
        .remove(&leaf_name)
        .ok_or_else(|| eyre::eyre!("Newick leaf '{leaf_name}' has no matching FASTA record"))?;
      Ok(Lock::new(Clade::new(Some(graph))))
    },
    (Some(l), Some(r)) => {
      let new_left = attach_graphs(&l, by_name)?;
      let new_right = attach_graphs(&r, by_name)?;
      Ok(Lock::new(Clade::from_children(None, &new_left, &new_right)))
    },
    _ => make_error!("Internal node with only one child encountered while attaching graphs to guide tree"),
  }
}

struct Parser<'a> {
  input: &'a str,
  pos: usize,
}

impl<'a> Parser<'a> {
  fn new(input: &'a str) -> Self {
    Self { input, pos: 0 }
  }

  fn eof(&self) -> bool {
    self.pos >= self.input.len()
  }

  fn remaining(&self) -> &str {
    self.input.get(self.pos..).unwrap_or("")
  }

  fn peek(&self) -> Option<char> {
    self.remaining().chars().next()
  }

  fn advance(&mut self) {
    if let Some(c) = self.peek() {
      self.pos += c.len_utf8();
    }
  }

  fn skip_ws(&mut self) {
    while let Some(c) = self.peek() {
      if c.is_whitespace() {
        self.advance();
      } else {
        break;
      }
    }
  }

  fn parse_subtree(&mut self) -> Result<Lock<Clade<Option<String>>>, Report> {
    self.skip_ws();
    if self.peek() == Some('(') {
      self.advance();
      self.skip_ws();
      let mut children = vec![self.parse_subtree()?];
      self.skip_ws();
      while self.peek() == Some(',') {
        self.advance();
        self.skip_ws();
        children.push(self.parse_subtree()?);
        self.skip_ws();
      }
      match self.peek() {
        Some(')') => self.advance(),
        Some(c) => {
          return make_error!("Newick: expected ')' or ',' at position {}, found '{}'", self.pos, c);
        },
        None => return make_error!("Newick: unexpected end of input, expected ')'"),
      }
      let _internal_label = self.parse_name();
      self.parse_branch_length()?;

      if children.len() != 2 {
        return make_error!(
          "Newick: internal node has {} children; only strictly bifurcating trees are supported",
          children.len()
        );
      }
      let mut iter = children.into_iter();
      let left = iter.next().unwrap();
      let right = iter.next().unwrap();
      Ok(Lock::new(Clade::from_children(None, &left, &right)))
    } else {
      let name = self
        .parse_name()
        .ok_or_else(|| eyre::eyre!("Newick: leaf without a name at position {}", self.pos))?;
      self.parse_branch_length()?;
      Ok(Lock::new(Clade::new(Some(name))))
    }
  }

  fn parse_name(&mut self) -> Option<String> {
    self.skip_ws();
    if self.peek() == Some('\'') {
      self.advance();
      let mut name = String::new();
      loop {
        match self.peek() {
          None => return None,
          Some('\'') => {
            self.advance();
            if self.peek() == Some('\'') {
              name.push('\'');
              self.advance();
            } else {
              break;
            }
          },
          Some(c) => {
            name.push(c);
            self.advance();
          },
        }
      }
      Some(name)
    } else {
      let mut name = String::new();
      while let Some(c) = self.peek() {
        if matches!(c, '(' | ')' | ',' | ':' | ';') || c.is_whitespace() {
          break;
        }
        name.push(c);
        self.advance();
      }
      if name.is_empty() { None } else { Some(name) }
    }
  }

  fn parse_branch_length(&mut self) -> Result<(), Report> {
    self.skip_ws();
    if self.peek() == Some(':') {
      self.advance();
      self.skip_ws();
      let start = self.pos;
      while let Some(c) = self.peek() {
        if c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '+' | '-') {
          self.advance();
        } else {
          break;
        }
      }
      if self.pos == start {
        return make_error!("Newick: expected a number after ':' at position {}", self.pos);
      }
    }
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::io::fasta::FastaRecord;
  use crate::pangraph::strand::Strand::Forward;
  use crate::representation::seq::Seq;
  use crate::tree::clade::WithName;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  impl WithName for Option<String> {
    fn name(&self) -> Option<&str> {
      self.as_deref()
    }
  }

  fn round_trip(input: &str) -> String {
    parse_newick(input).unwrap().read().to_newick()
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::basic_tree(              "((A,B),(C,D));",                   "((A,B),(C,D));")]
  #[case::branch_lengths(          "((A:0.1,B:0.2):0.3,C:0.4);",       "((A,B),C);"    )]
  #[case::internal_labels(         "((A,B)inner,C)root;",              "((A,B),C);"    )]
  #[case::whitespace_and_newlines( "(\n  (A , B) ,\n  ( C, D )\n);\n", "((A,B),(C,D));")]
  #[case::quoted_names(            "('foo bar',B);",                   "(foo bar,B);"  )]
  #[case::doubled_quote(           "('it''s',B);",                     "(it's,B);"     )]
  #[case::trailing_semicolon(      "((A,B),C)",                        "((A,B),C);"    )]
  #[case::scientific_notation(     "(A:1e-3,B:2.5E+2);",               "(A,B);"        )]
  #[trace]
  fn newick_round_trip(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(expected, round_trip(input));
  }

  #[rstest]
  #[case::empty("")]
  #[case::only_whitespace("   \n  ")]
  #[case::unbalanced_open("((A,B);")]
  #[case::unbalanced_close("A,B);")]
  #[case::multifurcation("(A,B,C);")]
  #[case::unifurcation("(A);")]
  #[case::leaf_no_name("(,B);")]
  #[case::trailing_garbage("(A,B);xyz")]
  #[case::missing_branch_number("(A:,B);")]
  fn newick_rejects_malformed_input(#[case] input: &str) {
    assert!(parse_newick(input).is_err(), "expected error for input: {input:?}");
  }

  fn singleton(name: &str, index: usize) -> Pangraph {
    Pangraph::singleton(
      FastaRecord {
        seq_name: name.to_owned(),
        desc: None,
        seq: Seq::from_str("ACGT"),
        index,
      },
      Forward,
      false,
    )
  }

  #[rstest]
  fn build_tree_from_newick_attaches_graphs() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("tree.nwk");
    std::fs::write(&path, "((A,B),C);").unwrap();

    let graphs = vec![singleton("A", 0), singleton("B", 1), singleton("C", 2)];
    let tree = build_tree_from_newick(&path, graphs).unwrap();

    assert!(tree.read().data.is_none());
    let leaves = collect_leaf_names(&tree);
    assert_eq!(leaves, vec!["A", "B", "C"]);
  }

  #[rstest]
  fn build_tree_from_newick_errors_on_extra_leaf_in_tree() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("tree.nwk");
    std::fs::write(&path, "((A,B),Z);").unwrap();

    let graphs = vec![singleton("A", 0), singleton("B", 1), singleton("C", 2)];
    let err = build_tree_from_newick(&path, graphs).unwrap_err().to_string();
    assert!(
      err.contains("'Z' has no matching FASTA record"),
      "unexpected error message: {err}"
    );
  }

  #[rstest]
  fn build_tree_from_newick_errors_on_missing_record_in_tree() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("tree.nwk");
    std::fs::write(&path, "(A,B);").unwrap();

    let graphs = vec![singleton("A", 0), singleton("B", 1), singleton("C", 2)];
    let err = build_tree_from_newick(&path, graphs).unwrap_err().to_string();
    assert!(err.contains("[C]"), "unexpected error message: {err}");
    assert!(
      err.contains("not present in the guide tree"),
      "unexpected error message: {err}"
    );
  }

  #[rstest]
  fn build_tree_from_newick_errors_on_duplicate_leaf() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("tree.nwk");
    std::fs::write(&path, "((A,B),A);").unwrap();

    let graphs = vec![singleton("A", 0), singleton("B", 1)];
    let err = build_tree_from_newick(&path, graphs).unwrap_err().to_string();
    assert!(
      err.contains("'A' has no matching FASTA record"),
      "unexpected error message: {err}"
    );
  }

  fn recurse_leaf_names(c: &Lock<Clade<Option<Pangraph>>>, out: &mut Vec<String>) {
    let g = c.read();
    match (&g.left, &g.right) {
      (None, None) => {
        let p = g.data.as_ref().unwrap();
        let name = p.paths.values().next().unwrap().name.clone().unwrap();
        out.push(name);
      },
      (Some(l), Some(r)) => {
        recurse_leaf_names(l, out);
        recurse_leaf_names(r, out);
      },
      _ => panic!("non-bifurcating internal node"),
    }
  }

  fn collect_leaf_names(tree: &Lock<Clade<Option<Pangraph>>>) -> Vec<String> {
    let mut out = vec![];
    recurse_leaf_names(tree, &mut out);
    out
  }
}
