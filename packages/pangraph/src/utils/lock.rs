use parking_lot::RawRwLock;
use parking_lot::lock_api::{ArcRwLockReadGuard, ArcRwLockWriteGuard, RwLock};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::hash::{Hash, Hasher};
use std::ops::Deref;
use std::sync::Arc;

#[derive(Debug)]
pub struct Lock<T>(Arc<RwLock<RawRwLock, T>>);

impl<T> Clone for Lock<T> {
  #[must_use]
  fn clone(&self) -> Self {
    Lock(Arc::clone(&self.0))
  }
}

impl<T> Lock<T> {
  pub fn new(x: T) -> Lock<T> {
    Lock(Arc::new(RwLock::new(x)))
  }

  pub fn read(&self) -> ArcRwLockReadGuard<RawRwLock, T> {
    self.0.read_arc()
  }

  pub fn write(&self) -> ArcRwLockWriteGuard<RawRwLock, T> {
    self.0.write_arc()
  }
}

impl<T> Deref for Lock<T> {
  type Target = Arc<RwLock<RawRwLock, T>>;

  fn deref(&self) -> &Self::Target {
    &self.0
  }
}

impl<T> PartialEq for Lock<T> {
  fn eq(&self, other: &Self) -> bool {
    Arc::ptr_eq(&self.0, &other.0)
  }
}

impl<T> Eq for Lock<T> {}

impl<T> Hash for Lock<T> {
  fn hash<H: Hasher>(&self, state: &mut H) {
    Arc::as_ptr(&self.0).hash(state);
  }
}

impl<T> Serialize for Lock<T>
where
  T: Serialize,
{
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: Serializer,
  {
    self.read().serialize(serializer)
  }
}

impl<'de, T> Deserialize<'de> for Lock<T>
where
  T: Deserialize<'de>,
{
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: Deserializer<'de>,
  {
    Deserialize::deserialize(deserializer).map(Lock::new)
  }
}
