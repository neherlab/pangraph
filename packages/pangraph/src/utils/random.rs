use rand::{Rng, SeedableRng, rng};
use rand_isaac::Isaac64Rng;

pub fn get_random_number_generator(seed: &Option<u64>) -> impl Rng + use<> {
  match seed {
    None => Isaac64Rng::from_rng(&mut rng()),
    Some(seed) => Isaac64Rng::seed_from_u64(*seed),
  }
}
