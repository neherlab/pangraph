#![allow(unsafe_code)]

mod buf;
mod index;
mod map;
mod options;
mod options_args;
mod ptr;

pub use index::Minimap2Index;
pub use map::{Minimap2Mapper, Minimap2PafRow, Minimap2PafRowSeq, Minimap2Reg, Minimap2RegExtra, Minimap2Result};
pub use options::Minimap2Preset;
pub use options_args::Minimap2Args;

// // SAFETY: index should be safe to pass between threads
// unsafe impl Send for Minimap2Index {}
