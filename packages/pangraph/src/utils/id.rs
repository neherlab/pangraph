use std::hash::{Hash, Hasher};
use twox_hash::XxHash64;

pub trait Id<T: From<usize>>: Hash {
  fn id(&self) -> T {
    id(self)
  }
}

pub fn id<In, Out>(x: In) -> Out
where
  In: Hash,
  Out: From<usize>,
{
  let mut hasher = XxHash64::with_seed(0);
  x.hash(&mut hasher);
  Out::from(hasher.finish() as usize)
}
