#pragma once

#include <mpfr.h>

#include <compare>

#include "Core/Config.h"

class HighPrecisionNumber {
 public:
  using DataType = mpfr_t;
  using PrecisionType = mpfr_prec_t;

 private:
  DataType value_;
  PrecisionType precision_;
  bool isInitialized_ = false;
  Real approxValue_ = 0.0;

 public:
  HighPrecisionNumber() : isInitialized_(false) {}
  explicit HighPrecisionNumber(const Real& val,
                               PrecisionType precision = defaultHighPrecision)
      : precision_(precision), isInitialized_(true) {
    mpfr_init2(value_, precision);
    mpfr_set_d(value_, val, MPFR_RNDN);
    approxValue_ = getDouble();
  }
  template <typename T>
    requires std::convertible_to<T, Real>
  explicit HighPrecisionNumber(const T& val,
                               PrecisionType precision = defaultHighPrecision)
      : HighPrecisionNumber(Real(val), precision) {}
  HighPrecisionNumber(const HighPrecisionNumber& rhs)
      : isInitialized_(rhs.isInitialized_) {
    rhs.checkInitialized();
    if (isInitialized_) {
      precision_ = rhs.precision_;
      mpfr_init2(value_, precision_);
      mpfr_set(value_, rhs.value_, MPFR_RNDN);
    }
    approxValue_ = getDouble();
  }
  const HighPrecisionNumber& operator=(const HighPrecisionNumber& rhs) {
    rhs.checkInitialized();
    if (this != &rhs) {
      if (isInitialized_ && rhs.isInitialized_) {
        mpfr_set(value_, rhs.value_, MPFR_RNDN);
      } else {
        isInitialized_ = true;
        precision_ = rhs.precision_;
        mpfr_init2(value_, precision_);
        mpfr_set(value_, rhs.value_, MPFR_RNDN);
      }
    }
    approxValue_ = getDouble();
    return *this;
  }

  template <typename T>
    requires std::convertible_to<T, Real>
  const HighPrecisionNumber& operator=(const T& rhs) {
    if (this != &rhs) {
      HighPrecisionNumber tmp(rhs);
      std::swap(*this, tmp);
    }
    approxValue_ = getDouble();
    return *this;
  }

  ~HighPrecisionNumber() { mpfr_clear(value_); }

 public:
  auto operator+(const HighPrecisionNumber& rhs) const -> HighPrecisionNumber;
  auto operator-(const HighPrecisionNumber& rhs) const -> HighPrecisionNumber;
  auto operator*(const HighPrecisionNumber& rhs) const -> HighPrecisionNumber;
  auto operator/(const HighPrecisionNumber& rhs) const -> HighPrecisionNumber;
  auto operator-() const -> HighPrecisionNumber;
  auto operator+=(const HighPrecisionNumber& rhs) -> HighPrecisionNumber&;
  auto operator/=(const HighPrecisionNumber& rhs) -> HighPrecisionNumber&;
  auto operator<=>(const HighPrecisionNumber& rhs) const
      -> std::strong_ordering;
  auto operator==(const HighPrecisionNumber& rhs) const -> bool;
  auto getDouble() const -> double;
  auto initialized() const -> bool;
  void print() const;
  friend HighPrecisionNumber sqrt(const HighPrecisionNumber& rhs);
  friend HighPrecisionNumber fabs(const HighPrecisionNumber& rhs);
  friend HighPrecisionNumber pow(const HighPrecisionNumber& base,
                                 const HighPrecisionNumber& exp);

 private:
  void checkInitialized() const {
    if (!isInitialized_) {
      throw std::runtime_error("HighPrecisionNumber is not initialized.");
    }
  }
  void checkInitialized(const HighPrecisionNumber& rhs) const {
    if (!isInitialized_ || !rhs.isInitialized_) {
      throw std::runtime_error("HighPrecisionNumber is not initialized");
    }
  }
};
