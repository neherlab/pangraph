use std::fmt::Display;

#[macro_export]
macro_rules! o {
  ($x:expr $(,)?) => {
    ToOwned::to_owned($x)
  };
}

pub fn quote(x: impl Display) -> String {
  format!("\"{x}\"")
}

pub fn quote_single(x: impl Display) -> String {
  format!("'{x}'")
}

pub fn str_slice_safe(s: &str, start: usize, end: usize) -> &str {
  if start >= s.len() || end == 0 || start > end {
    return "";
  }

  let safe_end = end.min(s.len());

  // Find valid UTF-8 boundaries
  let safe_start = s.chars().take(start).map(|c| c.len_utf8()).sum();
  let safe_end = s.chars().take(safe_end).map(|c| c.len_utf8()).sum();

  #[allow(clippy::string_slice)]
  &s[safe_start..safe_end]
}
