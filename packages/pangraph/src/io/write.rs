/// Adapt `std::fmt::Write` to `std::io::Write`
pub struct WriteAdapterFmtToIo<W: std::io::Write>(pub W);

impl<W: std::io::Write> std::fmt::Write for WriteAdapterFmtToIo<W> {
  fn write_str(&mut self, s: &str) -> std::fmt::Result {
    #[allow(clippy::map_err_ignore)]
    self.0.write_all(s.as_bytes()).map_err(|_| std::fmt::Error)?;
    Ok(())
  }
}

/// Adapt `std::io::Write` to `std::fmt::Write`
pub struct WriteAdapterIoToFmt<W: std::io::Write>(pub W);

impl<W: std::io::Write> std::fmt::Write for WriteAdapterIoToFmt<W> {
  fn write_str(&mut self, s: &str) -> std::fmt::Result {
    self.0.write_all(s.as_bytes()).map_err(|_| std::fmt::Error)
  }
}
