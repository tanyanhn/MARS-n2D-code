#ifndef MULTIVARIATEPOLYNOMIAL_H
#define MULTIVARIATEPOLYNOMIAL_H

#include <vector>

#include "Core/Config.h"
#include "Core/numlib.h"

/// Represent a degree Ndeg polynomial of Nvar variables.
/**
 * @See https://github.com/uekstrom/polymul/blob/master/polymul.h
 */
template <int Nvar, int Ndeg, class CoefType = Real>
class MultivariatePolynomial;

// this helper takes care of the reduction of order when finding derivatives
template <int Ndeg, int Var>
struct _diff_wrapper_2D;

// this helper takes care of the evaluation of a multivariate polynomial
template <int Ndeg, class CoefType>
struct _eval_wrapper_2D;

template <class CoefType, int Ndeg1, int Ndeg2>
struct _polynomial_multiplier_2D;

///  Template specialization
/**
 * Represent a bivariate polynomial
 */
template <int Ndeg, class CoefType>
class MultivariatePolynomial<2, Ndeg, CoefType> {
 public:
  template <class CT>
  using vector = std::vector<CT>;

  enum { Dim = 2, Order = Ndeg + 1 };
  enum { nBasis = binom(Ndeg + Dim, Dim) };

 public:
  /**
   * Construct a polynomial with constant coefficient c0
   *  and other terms zero.
   */
  MultivariatePolynomial(const CoefType &c0 = CoefType()) {
    for (int i = 0; i <= Ndeg; ++i)
      for (int j = 0; j <= Ndeg - i; ++j) coefs[i][j] = CoefType();
    coefs[0][0] = c0;
  }

  /**
   * Construct a polynomial according to the given coefficients
   *  ordered by graded lexicographical, for example 1 x y x^2 xy y^2 ..
   */
  MultivariatePolynomial(const std::initializer_list<CoefType> l) {
    auto q = l.begin();
    for (int deg = 0; deg <= Ndeg; ++deg)
      for (int i = deg; i >= 0; --i) coefs[i][deg - i] = *q++;
  }

 public:
  static constexpr int getOrder() { return Ndeg; }
  /**
   * return coefficient c[i][j]
   */
  const CoefType &coefficient(int i, int j) const { return coefs[i][j]; }

  CoefType &coefficient(int i, int j) { return coefs[i][j]; }

  const vector<vector<CoefType>> &getCoefs() const { return coefs; }

  vector<vector<CoefType>> &getCoefs() { return coefs; }

 public:
  /**
   * Evaluate the polynomial at a point p.
   */
  template <class T>
  CoefType operator()(const Vec<T, Dim> &p) const {
    return _eval_wrapper_2D<Ndeg, CoefType>::calc(*this, p);
  }

  /**
   * Calculate the derivative of this with respect to variable var.
   */
  template <int Var>
  auto diff() const {
    return _diff_wrapper_2D<Ndeg, Var>::calc(*this);
  }

 public:
  /**
   * Arithmetic
   */
  MultivariatePolynomial<Dim, Ndeg, CoefType> operator+(
      const MultivariatePolynomial<Dim, Ndeg, CoefType> &rhs) const {
    MultivariatePolynomial<Dim, Ndeg, CoefType> res;
    for (int i = 0; i <= Ndeg; ++i)
      for (int j = 0; j <= Ndeg - i; ++j)
        res.coefficient(i, j) = coefs[i][j] + rhs.coefficient(i, j);
    return res;
  }

  MultivariatePolynomial<Dim, Ndeg, CoefType> operator-() const {
    MultivariatePolynomial<Dim, Ndeg, CoefType> res;
    for (int i = 0; i <= Ndeg; ++i)
      for (int j = 0; j <= Ndeg - i; ++j) res.coefficient(i, j) = -coefs[i][j];
    return res;
  }

  MultivariatePolynomial<Dim, Ndeg, CoefType> operator-(
      const MultivariatePolynomial<Dim, Ndeg, CoefType> &rhs) const {
    MultivariatePolynomial<Dim, Ndeg, CoefType> res;
    for (int i = 0; i <= Ndeg; ++i)
      for (int j = 0; j <= Ndeg - i; ++j)
        res.coefficient(i, j) = coefs[i][j] - rhs.coefficient(i, j);
    return res;
  }

  template <int Ndeg1, int Ndeg2, class CT>
  friend MultivariatePolynomial<Dim, Ndeg1 + Ndeg2, CT> operator*(
      const MultivariatePolynomial<Dim, Ndeg1, CT> &lhs,
      const MultivariatePolynomial<Dim, Ndeg2, CT> &rhs);

  MultivariatePolynomial<Dim, Ndeg, CoefType> operator*(Real scalar) const {
    MultivariatePolynomial<Dim, Ndeg, CoefType> res;
    for (int i = 0; i <= Ndeg; ++i)
      for (int j = 0; j <= Ndeg - i; ++j)
        res.coefficient(i, j) = coefs[i][j] * scalar;
    return res;
  }

 public:
  /**
   * Output
   */
  friend std::ostream &operator<<(
      std::ostream &os, const MultivariatePolynomial<Dim, Ndeg, CoefType> &mp) {
    os << "A bivariate polynomial of degree " << Ndeg
       << ", with the coefficients :\n";
    for (int i = 0; i <= Ndeg; ++i) {
      for (int j = 0; j <= Ndeg - i; ++j) os << mp.coefs[i][j] << ", ";
      for (int j = Ndeg - i + 1; j <= Ndeg; ++j) os << "x, ";
      os << "\n";
    }
    return os;
  }

 protected:
  /**
   * A bivariate polynomial has the form
   *   p(x) = \sum_{i+j<=Ndeg} c[i][j]*x^iy^j.
   */
  vector<vector<CoefType>> coefs =
      vector<vector<CoefType>>(Order, vector<CoefType>(Order));
};

template <int Ndeg, class CoefType>
struct _eval_wrapper_2D {
  template <class T>
  static CoefType calc(const MultivariatePolynomial<2, Ndeg, CoefType> &mp,
                       const Vec<T, 2> &p) {
    // use Horner's scheme over coef[i]
    CoefType res = static_cast<CoefType>(0.0);
    T basis = static_cast<T>(1.0);
    for (int i = 0; i <= Ndeg; ++i) {
      CoefType res_i = static_cast<CoefType>(0.0);
      for (int j = Ndeg - i; j >= 0; --j) {
        res_i = res_i * p[1] + mp.coefficient(i, j);
      }
      res += res_i * basis;
      basis *= p[0];
    }
    return res;
  }
};

template <int Ndeg, int Var>
struct _diff_wrapper_2D {
  template <class CoefType>
  static MultivariatePolynomial<2, Ndeg - 1, CoefType> calc(
      const MultivariatePolynomial<2, Ndeg, CoefType> &mp) {
    MultivariatePolynomial<2, Ndeg - 1, CoefType> res;
    for (int i = 0; i <= Ndeg - 1; ++i)
      for (int j = 0; j <= Ndeg - 1 - i; ++j)
        /*        if constexpr(Var == 1)
                  res.coefficient(i,j) = (i+1) * mp.coefficient(i+1, j);
                else
                  res.coefficient(i,j) = (j+1) * mp.coefficient(i, j+1);*/
        res.coefficient(i, j) = (Var == 0) ? (i + 1) * mp.coefficient(i + 1, j)
                                           : (j + 1) * mp.coefficient(i, j + 1);
    return res;
  }
};

template <int Var>
struct _diff_wrapper_2D<0, Var> {
  template <class CoefType>
  static MultivariatePolynomial<2, 0, CoefType> calc(
      const MultivariatePolynomial<2, 0, CoefType> &mp0) {
    return MultivariatePolynomial<2, 0, CoefType>(0);
  }
};

template <int Ndeg1, int Ndeg2, class CoefType>
inline MultivariatePolynomial<2, Ndeg1 + Ndeg2, CoefType> operator*(
    const MultivariatePolynomial<2, Ndeg1, CoefType> &lhs,
    const MultivariatePolynomial<2, Ndeg2, CoefType> &rhs) {
  MultivariatePolynomial<2, Ndeg1 + Ndeg2, CoefType> res;
  _polynomial_multiplier_2D<CoefType, Ndeg1, Ndeg2>::mul(
      res.getCoefs(), lhs.getCoefs(), rhs.getCoefs());
  return res;
}

// Recursive template classes for multiplication.

template <class CoefType, int Ndeg1, int Ndeg2>
struct _polynomial_multiplier_2D {
  using VVC = std::vector<std::vector<CoefType>>;
  // _add_ the product between p1 and p2 to dst.
  static void mul(VVC &dst, const VVC &p1, const VVC &p2) {
    _polynomial_multiplier_2D<CoefType, Ndeg1, Ndeg2 - 1>::mul(dst, p1, p2);
    _polynomial_multiplier_2D<CoefType, Ndeg1, Ndeg2>::mul_fixOrder(dst, p1,
                                                                    p2);
  }
  // m2 is monomials in Nvar variables and of order Ndeg2.
  // _add_ the product of p1 and m2 to dst.
  static void mul_fixOrder(VVC &dst, const VVC &p1, const VVC &m2) {
    for (int i = 0; i <= Ndeg1; ++i)
      for (int j = 0; j <= Ndeg2; ++j)
        dst[Ndeg1 + Ndeg2 - i - j][i + j] +=
            p1[Ndeg1 - i][i] * m2[Ndeg2 - j][j];
    _polynomial_multiplier_2D<CoefType, Ndeg1 - 1, Ndeg2>::mul_fixOrder(dst, p1,
                                                                        m2);
  }
};

template <class CoefType, int Ndeg1>
struct _polynomial_multiplier_2D<CoefType, Ndeg1, 0> {
  using VVC = std::vector<std::vector<CoefType>>;
  static void mul(VVC &dst, const VVC &p1, const VVC &p2) {
    for (int i = 0; i <= Ndeg1; ++i)
      for (int j = 0; j <= Ndeg1 - i; ++j) dst[i][j] += p1[i][j] * p2[0][0];
  }
  static void mul_fixOrder(VVC &dst, const VVC &p1, const VVC &m2) {
    for (int i = 0; i <= Ndeg1; ++i)
      for (int j = 0; j <= Ndeg1 - i; ++j) dst[i][j] += p1[i][j] * m2[0][0];
  }
};

template <class CoefType, int Ndeg2>
struct _polynomial_multiplier_2D<CoefType, 0, Ndeg2> {
  using VVC = std::vector<std::vector<CoefType>>;
  static void mul(VVC &dst, const VVC &p1, const VVC &p2) {
    for (int i = 0; i < Ndeg2; ++i)
      for (int j = 0; j < Ndeg2 - i; ++j) dst[i][j] += p1[0][0] * p2[i][j];
  }
  static void mul_fixOrder(VVC &dst, const VVC &p1, const VVC &m2) {
    for (int j = 0; j <= Ndeg2; ++j)
      dst[Ndeg2 - j][j] += p1[0][0] * m2[Ndeg2 - j][j];
  }
};

template <class CoefType>
struct _polynomial_multiplier_2D<CoefType, 0, 0> {
  using VVC = std::vector<std::vector<CoefType>>;
  static void mul(VVC &dst, const VVC &p1, const VVC &p2) {
    dst[0][0] += p1[0][0] * p2[0][0];
  }
  static void mul_fixOrder(VVC &dst, const VVC &p1, const VVC &m2) {
    dst[0][0] += p1[0][0] * m2[0][0];
  }
};

#endif  // 　MULTIVARIATEPOLYNOMIAL_H
