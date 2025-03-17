
#include "HighPrecisionNumber.h"

auto HighPrecisionNumber::operator+(const HighPrecisionNumber& rhs) const
    -> HighPrecisionNumber {
  checkInitialized(rhs);
  HighPrecisionNumber result(0.0, std::max(precision_, rhs.precision_));
  mpfr_add(result.value_, value_, rhs.value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}

auto HighPrecisionNumber::operator+=(const HighPrecisionNumber& rhs)
    -> HighPrecisionNumber& {
  checkInitialized(rhs);
  auto result = *this + rhs;
  std::swap(*this, result);
  approxValue_ = getDouble();
  return *this;
}

auto HighPrecisionNumber::operator-(const HighPrecisionNumber& rhs) const
    -> HighPrecisionNumber {
  checkInitialized(rhs);
  HighPrecisionNumber result(0.0, std::max(precision_, rhs.precision_));
  mpfr_sub(result.value_, value_, rhs.value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}

auto HighPrecisionNumber::operator*(const HighPrecisionNumber& rhs) const
    -> HighPrecisionNumber {
  checkInitialized(rhs);
  HighPrecisionNumber result(0.0, std::max(precision_, rhs.precision_));
  mpfr_mul(result.value_, value_, rhs.value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}

auto HighPrecisionNumber::operator/(const HighPrecisionNumber& rhs) const
    -> HighPrecisionNumber {
  checkInitialized(rhs);
  HighPrecisionNumber result(0.0, std::max(precision_, rhs.precision_));
  mpfr_div(result.value_, value_, rhs.value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}

auto HighPrecisionNumber::operator/=(const HighPrecisionNumber& rhs)
    -> HighPrecisionNumber& {
  checkInitialized(rhs);
  auto result = *this / rhs;
  std::swap(*this, result);
  approxValue_ = getDouble();
  return *this;
}

auto HighPrecisionNumber::operator-() const -> HighPrecisionNumber {
  checkInitialized();
  HighPrecisionNumber result(0.0, precision_);
  mpfr_neg(result.value_, value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}

auto HighPrecisionNumber::operator<=>(const HighPrecisionNumber& rhs) const
    -> std::strong_ordering {
  checkInitialized(rhs);
  int cmp_result = mpfr_cmp(value_, rhs.value_);
  if (cmp_result == 0) {
    return std::strong_ordering::equal;
  }
  if (cmp_result < 0) {
    return std::strong_ordering::less;
  }
  return std::strong_ordering::greater;
}

auto HighPrecisionNumber::operator==(const HighPrecisionNumber& rhs) const
    -> bool {
  checkInitialized(rhs);
  return mpfr_equal_p(value_, rhs.value_);
}

auto HighPrecisionNumber::getDouble() const -> double {
  checkInitialized();
  return mpfr_get_d(value_, MPFR_RNDN);
}

void HighPrecisionNumber::print() const {
  checkInitialized();
  mpfr_printf("value = %Rf\n", value_);
}

HighPrecisionNumber sqrt(const HighPrecisionNumber& rhs) {
  rhs.checkInitialized();
  HighPrecisionNumber result(0.0, rhs.precision_);
  mpfr_sqrt(result.value_, rhs.value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}

HighPrecisionNumber fabs(const HighPrecisionNumber& rhs) {
  rhs.checkInitialized();
  HighPrecisionNumber result(0.0, rhs.precision_);
  mpfr_abs(result.value_, rhs.value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}

HighPrecisionNumber pow(const HighPrecisionNumber& base,
                        const HighPrecisionNumber& exp) {
  base.checkInitialized(exp);
  HighPrecisionNumber result(0.0, base.precision_);
  mpfr_pow(result.value_, base.value_, exp.value_, MPFR_RNDN);
  result.approxValue_ = result.getDouble();
  return result;
}
