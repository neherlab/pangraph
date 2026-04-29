#[macro_export]
macro_rules! pretty_assert_eq {
  ($left:expr, $right:expr) => {{
    pretty_assertions::assert_eq!(
      format!("{:#?}", $left).replace("\n", "\u{0085}"),
      format!("{:#?}", $right).replace("\n", "\u{0085}")
    );
  }};
}

#[macro_export]
macro_rules! assert_error {
  ($result:expr, $expected_message:expr) => {{
    let Err(error) = $result else {
      panic!("expected Err, got Ok");
    };
    let actual_message = $crate::utils::error::report_to_string(&error);
    pretty_assertions::assert_eq!($expected_message, actual_message);
  }};
}
