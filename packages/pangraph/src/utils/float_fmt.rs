use pretty_dtoa::{FmtFloatConfig, dtoa};
use std::sync::LazyLock;

static FLOAT_CONFIG: LazyLock<FmtFloatConfig> = LazyLock::new(|| {
  FmtFloatConfig::default()
    .force_no_e_notation()
    .add_point_zero(true)
    .max_significant_digits(3)
    .radix_point('.')
    .round()
});

pub fn float_to_significant_digits<F: Into<f64>>(weight: F, max_significant_digits: u8) -> String {
  dtoa(
    weight.into(),
    FLOAT_CONFIG.max_significant_digits(max_significant_digits),
  )
}
