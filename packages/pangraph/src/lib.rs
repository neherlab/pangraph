pub mod align;
pub mod circularize;
pub mod commands;
pub mod distance;
pub mod io;
pub mod pangraph;
pub mod reconsensus;
pub mod representation;
pub mod tree;
pub mod utils;

#[cfg(test)]
mod tests {
  use crate::utils::global_init::global_init;
  use ctor::ctor;

  #[ctor]
  fn init() {
    global_init();
  }
}
