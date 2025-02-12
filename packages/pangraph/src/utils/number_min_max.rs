use num_traits::Bounded;

pub trait MinMaxValueOfType: Bounded + Sized {
  /// Returns the minimum value for this integer type
  ///
  /// Example:
  ///
  /// ```rust
  /// use pangraph::utils::number_min_max::MinMaxValueOfType;
  /// use pretty_assertions::assert_eq;
  ///
  /// let x: i8 = 0;
  /// assert_eq!(x.min_value_of_type(), i8::MIN);
  ///
  /// let y: usize = 123456;
  /// assert_eq!(y.min_value_of_type(), usize::MIN);
  /// ```
  #[must_use]
  fn min_value_of_type(&self) -> Self {
    <Self as Bounded>::min_value()
  }

  /// Returns the maximum value for this integer type
  ///
  /// Example:
  ///
  /// ```rust
  /// use pangraph::utils::number_min_max::MinMaxValueOfType;
  /// use pretty_assertions::assert_eq;
  ///
  /// let x: i8 = 0;
  /// assert_eq!(x.max_value_of_type(), i8::MAX);
  ///
  /// let y: usize = 123456;
  /// assert_eq!(y.max_value_of_type(), usize::MAX);
  /// ```
  #[must_use]
  fn max_value_of_type(&self) -> Self {
    <Self as Bounded>::max_value()
  }
}

pub trait IsMinMaxValueOfType: MinMaxValueOfType + Eq {
  /// Is the current value the minimum value for this integer type?
  ///
  /// Example:
  ///
  /// ```rust
  /// use pangraph::utils::number_min_max::IsMinMaxValueOfType;
  ///
  /// let x: i8 = i8::MIN;
  /// assert!(x.is_min_value_of_type());
  ///
  /// let y: i8 = 5;
  /// assert!(!y.is_min_value_of_type());
  /// ```
  #[must_use]
  fn is_min_value_of_type(&self) -> bool {
    self == &self.min_value_of_type()
  }

  /// Is the current value the maximum value for this integer type?
  ///
  /// Example:
  ///
  /// ```rust
  /// use pangraph::utils::number_min_max::IsMinMaxValueOfType;
  ///
  /// let x: i8 = i8::MAX;
  /// assert!(x.is_max_value_of_type());
  ///
  /// let y: i8 = 5;
  /// assert!(!y.is_max_value_of_type());
  /// ```
  #[must_use]
  fn is_max_value_of_type(&self) -> bool {
    self == &self.max_value_of_type()
  }
}

impl<T: Bounded> MinMaxValueOfType for T {}

impl<T: MinMaxValueOfType + Eq> IsMinMaxValueOfType for T {}
