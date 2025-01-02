#[cfg(test)]
mod tests {
  use ctor::ctor;
  use pangraph::utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
