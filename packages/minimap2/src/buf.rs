#![allow(unsafe_code)]
use eyre::{eyre, Report};
use minimap2_sys::{mm_tbuf_destroy, mm_tbuf_init, mm_tbuf_t};
use std::ptr::NonNull;

/// NOTE: this buffer is not thread-safe. Each thread must use a dedicated buffer.
#[must_use]
#[derive(Debug)]
pub struct Minimap2Buffer {
  buffer: NonNull<mm_tbuf_t>,
}

impl Minimap2Buffer {
  pub fn new() -> Result<Self, Report> {
    // SAFETY: calling unsafe function. Need to ensure buffer is not null.
    unsafe {
      NonNull::new(mm_tbuf_init())
        .map(|buffer| Self { buffer })
        .ok_or_else(|| eyre!("Buffer allocation failed"))
    }
  }

  pub fn get(&self) -> *const mm_tbuf_t {
    self.buffer.as_ptr()
  }

  pub fn get_mut(&self) -> *mut mm_tbuf_t {
    self.buffer.as_ptr()
  }
}

impl Drop for Minimap2Buffer {
  fn drop(&mut self) {
    // SAFETY: calling unsafe function. Need to ensure buffer is not null.
    unsafe {
      mm_tbuf_destroy(self.buffer.as_ptr());
    }
  }
}
