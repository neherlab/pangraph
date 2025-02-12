use std::ffi::{CString, NulError};
use std::mem::MaybeUninit;
use std::os::raw::c_char;

pub fn to_c_char_ptr(input: &str) -> Result<*const c_char, NulError> {
  CString::new(input).map(|cs| cs.into_raw().cast_const())
}

pub fn to_c_char_ptrs(input: &[&str]) -> Result<*mut *const c_char, NulError> {
  let c_strings: Vec<CString> = input.iter().map(|&s| CString::new(s)).collect::<Result<_, _>>()?;
  let c_ptrs: Vec<*const c_char> = c_strings.iter().map(|s| s.as_ptr()).collect();
  Ok(c_ptrs.into_boxed_slice().as_mut_ptr())
}

/// Check if the pointer associated with a `MaybeUninit<T>` is null.
pub fn is_null<T>(maybe_uninit: &MaybeUninit<T>) -> bool {
  maybe_uninit.as_ptr().is_null()
}
