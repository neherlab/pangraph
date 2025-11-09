use crate::io::compression::Decompressor;
use crate::io::concat::Concat;
use crate::io::file::{create_file_or_stdout, open_file_or_stdin, open_stdin};
use crate::make_error;
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use crate::utils::string::{quote_single, str_slice_safe};
use color_eyre::{Section, SectionExt};
use eyre::{Context, Report};
use itertools::Itertools;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

const ALPHABET: &[char] = &[
  'A', 'C', 'G', 'T', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'N', '-',
];

#[derive(Clone, Default, Debug, Deserialize, Serialize, Eq, PartialEq)]
#[serde(rename_all = "camelCase")]
pub struct FastaRecord {
  pub seq_name: String,
  pub desc: Option<String>,
  pub seq: Seq,
  pub index: usize,
}

impl FastaRecord {
  pub fn new() -> Self {
    Self::default()
  }

  pub fn clear(&mut self) {
    self.seq_name.clear();
    self.desc = None;
    self.seq.clear();
    self.index = 0;
  }

  pub fn is_empty(&self) -> bool {
    self.seq_name.is_empty() && self.seq_name.is_empty() && self.desc.is_none() && self.index == 0
  }

  pub fn header(&self) -> String {
    match &self.desc {
      Some(desc) => format!(">{} {}", self.seq_name, desc),
      None => format!(">{}", self.seq_name),
    }
  }
}

pub struct FastaReader<'a> {
  reader: Box<dyn BufRead + 'a>,
  line: String,
  n_lines: usize,
  n_chars: usize,
  index: usize,
  disallowed_chars: Vec<char>,
}

impl<'a> FastaReader<'a> {
  pub fn new(reader: Box<dyn BufRead + 'a>) -> Self {
    Self {
      reader,
      line: String::new(),
      n_lines: 0,
      n_chars: 0,
      index: 0,
      disallowed_chars: Vec::new(),
    }
  }

  #[must_use]
  pub fn with_disallowed_chars(mut self, disallowed_chars: Vec<char>) -> Self {
    self.disallowed_chars = disallowed_chars;
    self
  }

  pub fn from_str(contents: &'a impl AsRef<str>) -> Result<Self, Report> {
    let reader = contents.as_ref().as_bytes();
    Ok(Self::new(Box::new(reader)))
  }

  pub fn from_str_and_path(contents: &'static str, filepath: impl AsRef<Path>) -> Result<Self, Report> {
    let decompressor = Decompressor::from_str_and_path(contents, filepath)?;
    let reader = BufReader::new(decompressor);
    Ok(Self::new(Box::new(reader)))
  }

  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Self::from_paths(&[filepath])
  }

  /// Reads multiple files sequentially given a set of paths
  pub fn from_paths<P: AsRef<Path>>(filepaths: &[P]) -> Result<Self, Report> {
    if filepaths.is_empty() {
      info!("Reading input fasta from standard input");
      return Ok(Self::new(open_stdin()?));
    }

    let readers: Vec<Box<dyn BufRead + 'a>> = filepaths
      .iter()
      .map(|filepath| -> Result<Box<dyn BufRead + 'a>, Report> { open_file_or_stdin(&Some(filepath)) })
      .collect::<Result<Vec<Box<dyn BufRead + 'a>>, Report>>()?;

    let concat = Concat::with_delimiter(readers.into_iter(), Some(b"\n".to_vec()));
    let concat_buf = BufReader::new(concat);

    Ok(Self::new(Box::new(concat_buf)))
  }

  #[allow(clippy::string_slice)]
  pub fn read(&mut self, record: &mut FastaRecord) -> Result<(), Report> {
    record.clear();

    if self.line.is_empty() {
      // Read lines until we find the next record or EOF
      loop {
        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
          if self.index > 0 {
            // We have read at least one record  by the end of input - this is normal operation mode
            return Ok(());
          }

          if self.index == 0 && self.n_chars == 0 {
            // We have read no records and no non-whitespace characters by the end of input
            warn!(
              "FASTA input is empty or consists entirely from whitespace: this is allowed but might not be what's intended"
            );
            return Ok(());
          }

          // We have read some characters, but no records detected by the end of input
          return make_error!(
            "FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found"
          );
        }

        let trimmed = self.line.trim();
        self.n_lines += 1;
        self.n_chars += trimmed.len();
        if trimmed.starts_with('>') {
          break; // Found the header of the next record
        }
      }
    }

    let header_line = self.line.trim();
    let (name, desc) = header_line[1..]
      .split_once(' ')
      .unwrap_or_else(|| (&header_line[1..], ""));
    record.seq_name = name.to_owned();
    record.desc = if desc.is_empty() { None } else { Some(desc.to_owned()) };
    record.index = self.index;
    self.index += 1;

    // Read sequence lines until the next header or EOF
    self.line.clear();
    while self.reader.read_line(&mut self.line)? > 0 {
      let trimmed = self.line.trim();
      self.n_lines += 1;
      self.n_chars += trimmed.len();
      if trimmed.starts_with('>') {
        // We have reached the next record
        break;
      }

      record.seq.reserve(trimmed.len());
      for (i, c) in trimmed.chars().enumerate() {
        let c_upper = c.to_ascii_uppercase();
        if !ALPHABET.contains(&c_upper) {
          return make_error!("FASTA input is incorrect: character \"{c}\" is not in the alphabet",)
            .with_section(|| {
              str_slice_safe(trimmed, i.saturating_sub(10), i.saturating_add(10))
                .to_owned()
                .header("Sequence fragment:")
            })
            .with_section(|| {
              ALPHABET
                .iter()
                .map(quote_single)
                .join(", ")
                .header("Expected characters: ")
            })
            .wrap_err_with(|| format!("When processing sequence #{}: \"{}\"", self.index, record.header()));
        }
        if self.disallowed_chars.contains(&c_upper) {
          return make_error!("FASTA input is incorrect: character \"{c}\" is disallowed",)
            .with_section(|| {
              str_slice_safe(trimmed, i.saturating_sub(10), i.saturating_add(10))
                .to_owned()
                .header("Sequence fragment:")
            })
            .with_section(|| {
              self
                .disallowed_chars
                .iter()
                .map(quote_single)
                .join(", ")
                .header("Disallowed characters: ")
            })
            .wrap_err_with(|| format!("When processing sequence #{}: \"{}\"", self.index, record.header()));
        }
        record.seq.push(AsciiChar::from(c_upper));
      }

      self.line.clear();
    }

    Ok(())
  }
}

pub fn read_one_fasta(filepath: impl AsRef<Path>) -> Result<FastaRecord, Report> {
  let filepath = filepath.as_ref();
  let mut reader = FastaReader::from_path(filepath)?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

pub fn read_many_fasta<P: AsRef<Path>>(filepaths: &[P]) -> Result<Vec<FastaRecord>, Report> {
  read_many_fasta_impl(filepaths, &[])
}

pub fn read_many_fasta_with_disallowed_chars<P: AsRef<Path>>(
  filepaths: &[P],
  disallowed_chars: &[char],
) -> Result<Vec<FastaRecord>, Report> {
  read_many_fasta_impl(filepaths, disallowed_chars)
}

fn read_many_fasta_impl<P: AsRef<Path>>(
  filepaths: &[P],
  disallowed_chars: &[char],
) -> Result<Vec<FastaRecord>, Report> {
  let mut reader = FastaReader::from_paths(filepaths)?.with_disallowed_chars(disallowed_chars.to_vec());
  let mut fasta_records = Vec::<FastaRecord>::new();

  loop {
    let mut record = FastaRecord::default();
    reader.read(&mut record)?;
    if record.is_empty() {
      break;
    }
    fasta_records.push(record);
  }

  Ok(fasta_records)
}

pub fn read_one_fasta_str(contents: impl AsRef<str>) -> Result<FastaRecord, Report> {
  let mut reader = FastaReader::from_str(&contents)?;
  let mut record = FastaRecord::default();
  reader.read(&mut record)?;
  Ok(record)
}

pub fn read_many_fasta_str(contents: impl AsRef<str>) -> Result<Vec<FastaRecord>, Report> {
  let mut reader = FastaReader::from_str(&contents)?;
  let mut fasta_records = Vec::<FastaRecord>::new();

  loop {
    let mut record = FastaRecord::default();
    reader.read(&mut record)?;
    if record.is_empty() {
      break;
    }
    fasta_records.push(record);
  }

  Ok(fasta_records)
}

// Writes sequences into given fasta file
pub struct FastaWriter {
  writer: Box<dyn Write>,
}

impl FastaWriter {
  pub fn new(writer: Box<dyn Write>) -> Self {
    Self { writer }
  }

  pub fn from_path(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Ok(Self::new(create_file_or_stdout(filepath)?))
  }

  pub fn write(&mut self, seq_name: impl AsRef<str>, desc: &Option<String>, seq: &Seq) -> Result<(), Report> {
    self.writer.write_all(b">")?;
    self.writer.write_all(seq_name.as_ref().as_bytes())?;

    if let Some(desc) = desc {
      self.writer.write_all(b" ")?;
      self.writer.write_all(desc.as_bytes())?;
    }

    self.writer.write_all(b"\n")?;
    self.writer.write_all(seq.as_ref())?;
    self.writer.write_all(b"\n")?;
    Ok(())
  }

  pub fn flush(&mut self) -> Result<(), Report> {
    self.writer.flush()?;
    Ok(())
  }
}

pub fn write_one_fasta(
  filepath: impl AsRef<Path>,
  seq_name: impl AsRef<str>,
  desc: &Option<String>,
  seq: &Seq,
) -> Result<(), Report> {
  let mut writer = FastaWriter::from_path(&filepath)?;
  writer.write(seq_name, desc, seq)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::utils::error::report_to_string;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::io::Cursor;

  #[test]
  fn test_fasta_reader_fail_on_non_fasta() {
    let data =
      b"This is not a valid FASTA string.\nIt is not empty, and not entirely whitespace\nbut does not contain 'greater than' character.\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));
    let mut record = FastaRecord::new();
    assert_eq!(
      reader.read(&mut record).unwrap_err().to_string(),
      "FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found"
    );
  }

  #[test]
  fn test_fasta_reader_fail_on_unknown_char() {
    let data = b">seq1\nACGT%ACGT\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));
    let mut record = FastaRecord::new();
    let actual = report_to_string(&reader.read(&mut record).unwrap_err());
    let expected =
      r#"When processing sequence #1: ">seq1": FASTA input is incorrect: character "%" is not in the alphabet"#;
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_fasta_reader_read_empty() {
    let data = b"";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert!(record.is_empty());
  }

  #[test]
  fn test_fasta_reader_read_whitespace_only() {
    let data = b"\n \n \n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert!(record.is_empty());
  }

  #[test]
  fn test_fasta_reader_read_single_record() {
    let data = b">seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_with_leading_newline() {
    let data = b"\n>seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_with_multiple_leading_newlines() {
    let data = b"\n\n\n>seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_without_trailing_newline() {
    let data = b">seq1\nATCG";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_multiple_records() {
    let data = b">seq1\nATCG\n>seq2\nGCTA\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record1 = FastaRecord::new();
    reader.read(&mut record1).unwrap();

    let mut record2 = FastaRecord::new();
    reader.read(&mut record2).unwrap();

    assert_eq!(record1.seq_name, "seq1");
    assert_eq!(record1.seq, "ATCG");
    assert_eq!(record1.index, 0);

    assert_eq!(record2.seq_name, "seq2");
    assert_eq!(record2.seq, "GCTA");
    assert_eq!(record2.index, 1);
  }

  #[test]
  fn test_fasta_reader_read_empty_lines_between_records() {
    let data = b"\n>seq1\n\nATCG\n\n\n>seq2\nGCTA\n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record1 = FastaRecord::new();
    reader.read(&mut record1).unwrap();

    let mut record2 = FastaRecord::new();
    reader.read(&mut record2).unwrap();

    assert_eq!(record1.seq_name, "seq1");
    assert_eq!(record1.seq, "ATCG");
    assert_eq!(record1.index, 0);

    assert_eq!(record2.seq_name, "seq2");
    assert_eq!(record2.seq, "GCTA");
    assert_eq!(record2.index, 1);
  }

  #[test]
  fn test_fasta_reader_read_with_trailing_newline() {
    let data = b">seq1\nATCG\n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_example_1() {
    let data = b"\n\n>a\nACGCTCGATC\n\n>b\nCCGCGC";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: "ACGCTCGATC".into(),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: "CCGCGC".into(),
        index: 1,
      }
    );
  }

  #[test]
  fn test_fasta_reader_example_2() {
    let data = b">a\nACGCTCGATC\n>b\nCCGCGC\n>c";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: "ACGCTCGATC".into(),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: "CCGCGC".into(),
        index: 1,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("c"),
        desc: None,
        seq: "".into(),
        index: 2,
      }
    );
  }

  #[test]
  fn test_fasta_reader_example_3() {
    let data = b">a\nACGCTCGATC\n>b\n>c\nCCGCGC";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)));

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: "ACGCTCGATC".into(),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: "".into(),
        index: 1,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("c"),
        desc: None,
        seq: "CCGCGC".into(),
        index: 2,
      }
    );
  }

  #[test]
  fn test_fasta_reader_name_desc() -> Result<(), Report> {
    let actual = read_many_fasta_str(indoc! {r#"
      >Identifier Description
      ACGT
      >Identifier Description with spaces
      ACGT


    "#})?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description")),
        seq: "ACGT".into(),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description with spaces")),
        seq: "ACGT".into(),
        index: 1,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_dedent_nuc() -> Result<(), Report> {
    let actual = read_many_fasta_str(indoc! {r#"
      >FluBuster-001
      ACAGCCATGTATTG--
      >CommonCold-AB
      ACATCCCTGTA-TG--
      >Ecoli/Joke/2024|XD
      ACATCGCCNNA--GAC

      >Sniffles-B
      GCATCCCTGTA-NG--
      >StrawberryYogurtCulture|🍓
      CCGGCCATGTATTG--
      > SneezeC-19
      CCGGCGATGTRTTG--
        >MisindentedVirus|D-skew
        TCGGCCGTGTRTTG--
    "#})?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("FluBuster-001"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("CommonCold-AB"),
        desc: None,
        seq: "ACATCCCTGTA-TG--".into(),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("Ecoli/Joke/2024|XD"),
        desc: None,
        seq: "ACATCGCCNNA--GAC".into(),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("Sniffles-B"),
        desc: None,
        seq: "GCATCCCTGTA-NG--".into(),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("StrawberryYogurtCulture|🍓"),
        desc: None,
        seq: "CCGGCCATGTATTG--".into(),
        index: 4,
      },
      FastaRecord {
        seq_name: o!(""),
        desc: Some(o!("SneezeC-19")),
        seq: "CCGGCGATGTRTTG--".into(),
        index: 5,
      },
      FastaRecord {
        seq_name: o!("MisindentedVirus|D-skew"),
        desc: None,
        seq: "TCGGCCGTGTRTTG--".into(),
        index: 6,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_multiline_and_skewed_indentation() -> Result<(), Report> {
    let actual = read_many_fasta_str(indoc! {r#"
      >MixedCaseSeq
      aCaGcCAtGtAtTG--
      >LowercaseSeq
      acagccatgtattg--
      >UppercaseSeq
      ACAGCCATGTATTG--
      >MultilineSeq
      ACAGCC
      ATGT
      ATTG--
      >SkewedIndentSeq
        ACAGCC
      ATGTATTG
       ATTG--
    "#})?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("MixedCaseSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("LowercaseSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("UppercaseSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("MultilineSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("SkewedIndentSeq"),
        desc: None,
        seq: "ACAGCCATGTATTGATTG--".into(),
        index: 4,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_disallowed_gap_characters() {
    let data = b">seq1\nACGT\n>seq2\nAC-GT\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data))).with_disallowed_chars(vec!['-']);

    // First record should succeed (no gaps)
    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();
    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ACGT");

    // Second record should fail (contains gap)
    let mut record2 = FastaRecord::new();
    let result = reader.read(&mut record2);
    assert!(result.is_err());
    
    let error_message = report_to_string(&result.unwrap_err());
    assert!(error_message.contains("character \"-\" is disallowed"));
  }
}
