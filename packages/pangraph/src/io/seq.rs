use crate::make_error;
use eyre::Report;
use rand::seq::SliceRandom;
use rand::Rng;

pub fn complement(nuc: char) -> Result<char, Report> {
  Ok(match nuc {
    'A' => 'T',
    'C' => 'G',
    'G' => 'C',
    'T' => 'A',
    'Y' => 'R',
    'R' => 'Y',
    'W' => 'W',
    'S' => 'S',
    'K' => 'M',
    'M' => 'K',
    'D' => 'H',
    'V' => 'B',
    'H' => 'D',
    'B' => 'V',
    'N' => 'N',
    '-' => '-',
    _ => return make_error!("Unknown nucleotide character: '{nuc}'"),
  })
}

pub fn reverse_complement(seq: impl AsRef<str>) -> Result<String, Report> {
  seq
    .as_ref()
    .chars()
    .rev()
    .map(complement)
    .collect::<Result<String, Report>>()
}

pub fn generate_random_nuc_sequence(length: usize, rng: &mut impl Rng) -> String {
  const CHOICES: [char; 4] = ['A', 'C', 'G', 'T'];
  (0..length)
    .map(|_| CHOICES.choose(rng).expect("choosing from an empty set"))
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
    assert_eq!(reverse_complement("ATCG")?, "CGAT");
    assert_eq!(reverse_complement("GATTACA")?, "TGTAATC");
    assert_eq!(reverse_complement("ACGTYRWSKM")?, "KMSWYRACGT");
    assert_eq!(reverse_complement("N-")?, "-N");
    Ok(())
  }

  #[test]
  fn test_reverse_complement_invalid() {
    assert_eq!(
      reverse_complement("ATCGX").unwrap_err().to_string(),
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
