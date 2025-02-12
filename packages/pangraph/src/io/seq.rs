use crate::make_error;
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use eyre::Report;
use rand::seq::SliceRandom;
use rand::Rng;

pub fn complement(nuc: &AsciiChar) -> Result<AsciiChar, Report> {
  Ok(match *nuc {
    AsciiChar(b'A') => AsciiChar(b'T'),
    AsciiChar(b'C') => AsciiChar(b'G'),
    AsciiChar(b'G') => AsciiChar(b'C'),
    AsciiChar(b'T') => AsciiChar(b'A'),
    AsciiChar(b'Y') => AsciiChar(b'R'),
    AsciiChar(b'R') => AsciiChar(b'Y'),
    AsciiChar(b'W') => AsciiChar(b'W'),
    AsciiChar(b'S') => AsciiChar(b'S'),
    AsciiChar(b'K') => AsciiChar(b'M'),
    AsciiChar(b'M') => AsciiChar(b'K'),
    AsciiChar(b'D') => AsciiChar(b'H'),
    AsciiChar(b'V') => AsciiChar(b'B'),
    AsciiChar(b'H') => AsciiChar(b'D'),
    AsciiChar(b'B') => AsciiChar(b'V'),
    AsciiChar(b'N') => AsciiChar(b'N'),
    AsciiChar(b'-') => AsciiChar(b'-'),
    _ => return make_error!("Unknown nucleotide character: '{nuc}'"),
  })
}

pub fn reverse_complement(seq: &Seq) -> Result<Seq, Report> {
  seq.iter().rev().map(complement).collect::<Result<Seq, Report>>()
}

pub fn generate_random_nuc_sequence(length: usize, rng: &mut impl Rng) -> Seq {
  const CHOICES: [char; 4] = ['A', 'C', 'G', 'T'];
  (0..length)
    .map(|_| CHOICES.choose(rng).expect("choosing from an empty set"))
    .map(|c| AsciiChar::from(*c))
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::utils::random::get_random_number_generator;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_reverse_complement() -> Result<(), Report> {
    assert_eq!(reverse_complement(&Seq::from("ATCG"))?, "CGAT");
    assert_eq!(reverse_complement(&Seq::from("GATTACA"))?, "TGTAATC");
    assert_eq!(reverse_complement(&Seq::from("ACGTYRWSKM"))?, "KMSWYRACGT");
    assert_eq!(reverse_complement(&Seq::from("N-"))?, "-N");
    Ok(())
  }

  #[test]
  fn test_reverse_complement_invalid() {
    assert_eq!(
      reverse_complement(&Seq::from("ATCGX")).unwrap_err().to_string(),
      "Unknown nucleotide character: 'X'"
    );
  }

  #[test]
  fn test_generate_random_nuc_sequence() -> Result<(), Report> {
    let mut rng = get_random_number_generator(&Some(0));
    assert_eq!(
      generate_random_nuc_sequence(123, &mut rng),
      "GGGGCGGACCAATCTCCCTACTGCCAGCGCTCCGGCCAATCGAGGCCCCCAAAATTTGGATGCCATGACCGCAATTATGAACATACTCCGCTGTGTAACCTTTGGTCAGAGCTCGCGGAGAAT"
    );
    Ok(())
  }
}
