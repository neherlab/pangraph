use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Edits {
  subs: Vec<Sub>,
  dels: Vec<Del>,
  inss: Vec<Ins>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Sub {
  pos: usize,
  alt: u8,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Del {
  pos: usize,
  len: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Ins {
  pos: usize,
  ins: String,
}
