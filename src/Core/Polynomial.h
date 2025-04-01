#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "numlib.h"
#include "static_for.h"

// this helper takes care of the reduction of order when finding derivatives
template <int Order>
struct der_wrapper;

///
/**
 */
template <int Order, class CoefType>
class Polynomial {
 protected:
  CoefType coefs[Order];

 public:
  enum { _Order = Order };

  Polynomial(const CoefType &_c = CoefType()) {
    coefs[0] = _c;
    for (int i = 1; i < Order; coefs[i++] = CoefType())
      ;
  }

  Polynomial(std::initializer_list<CoefType> l) {
    auto q = l.begin();
    for (int i = 0; i < Order; ++i) coefs[i] = *q++;
  }

  template <class CoefType2>
  explicit Polynomial(const Polynomial<Order, CoefType2> &p) {
    for (int i = 0; i < Order; ++i) coefs[i] = p[i];
  }

 public:
  const CoefType &operator[](int k) const { return coefs[k]; }

  CoefType &operator[](int k) { return coefs[k]; }

  ///
  /**
     Evaluate the polynomial at x.
   */
  template <class T2>
#ifdef OPTNONE
  __attribute__((optnone))
#endif  // OPTNONE
  CoefType
  operator()(const T2 &x) const {
    CoefType val = 0;
    for (int d = Order - 1; d >= 0; --d) {
      val = val * x + coefs[d];
    }
    return val;
  }

  auto der() const { return der_wrapper<Order>::calc(*this); }

  ///
  /**
     The constant term is always 0.
   */
  Polynomial<Order + 1, CoefType> prim() const {
    Polynomial<Order + 1, CoefType> res;
    for (int d = 1; d < Order + 1; ++d) res[d] = coefs[d - 1] / d;
    return res;
  }

  ///
  /**
     Re-express in the y = (x-x0) if reflect = false,
     in y = (x0-x) if reflect = true.
     input f(x) = f(x0 + y) or f(x0 - y) -> return g(y)
   */
  template <class T2>
  Polynomial<Order, CoefType> translate(const T2 &x0,
                                        bool reflect = false) const {
    // note : the Taylor expansion of a polynomial of order D
    // up to D-th order gives the polynomial itself
    Polynomial<Order, CoefType> res;
    res[0] = (*this)(x0);
    auto calc = [&x0, &res](auto _i, const auto &poly) {
      auto dp = poly.der();
      res[decltype(_i)::value] = dp(x0) / factorial(decltype(_i)::value);
      return dp;
    };
    static_for<1, Order>::execute(calc, *this);
    if (reflect) {
      for (int d = 1; d < Order; d += 2)  // negate the odd terms
        res.coefs[d] = -res.coefs[d];
    }
    return res;
  }

  // arithmetics
 public:
  Polynomial<Order, CoefType> operator+(
      const Polynomial<Order, CoefType> &rhs) const {
    Polynomial<Order, CoefType> res;
    for (int d = 0; d < Order; d++) res[d] = coefs[d] + rhs[d];
    return res;
  }

  Polynomial<Order, CoefType> operator-() const {
    Polynomial<Order, CoefType> res;
    for (int d = 0; d < Order; d++) res[d] = -coefs[d];
    return res;
  }

  Polynomial<Order, CoefType> operator-(
      const Polynomial<Order, CoefType> &rhs) const {
    Polynomial<Order, CoefType> res;
    for (int d = 0; d < Order; d++) res[d] = coefs[d] - rhs[d];
    return res;
  }

  template <int Ord1, int Ord2, class CT>
  friend Polynomial<Ord1 + Ord2 - 1, CT> operator*(
      const Polynomial<Ord1, CT> &lhs, const Polynomial<Ord2, CT> &rhs);

  Polynomial<Order, CoefType> operator*(Real scalar) const {
    Polynomial<Order, CoefType> res;
    for (int d = 0; d < Order; d++) res[d] = coefs[d] * scalar;
    return res;
  }

  Polynomial<Order, CoefType> operator/(Real scalar) const {
    Polynomial<Order, CoefType> res;
    for (int d = 0; d < Order; d++) res.coefs[d] = coefs[d] / scalar;
    return res;
  }

  // CoefType diffInClosedPoint(Real t0, Real t1) const {  }
};

template <class CoefType>
class Polynomial<0, CoefType>;  // 0 order is meaningless

//==================================================

template <int Order>
struct der_wrapper {
  template <class CoefType>
  static Polynomial<Order - 1, CoefType> calc(
      const Polynomial<Order, CoefType> &P) {
    Polynomial<Order - 1, CoefType> res;
    for (int d = 0; d < Order - 1; d++) res[d] = P[d + 1] * (d + 1);
    return res;
  }
};

template <>
struct der_wrapper<1> {
  template <class CoefType>
  static Polynomial<1, CoefType> calc(const Polynomial<1, CoefType> &P) {
    return Polynomial<1, CoefType>(0);
  }
};

template <int Order, class T>
struct ipow_wrapper<Polynomial<Order, T>, 0, 0> {
  static Polynomial<1, T> calc(const Polynomial<Order, T> &) {
    return Polynomial<1, T>(1);
  }
};

template <int Ord1, int Ord2, class CT>
inline Polynomial<Ord1 + Ord2 - 1, CT> operator*(
    const Polynomial<Ord1, CT> &lhs, const Polynomial<Ord2, CT> &rhs) {
  Polynomial<Ord1 + Ord2 - 1, CT> res;  // res is zeroized
  for (int k = 0; k < Ord1 + Ord2 - 1; ++k)
    for (int j = std::max(0, k - Ord2 + 1); j <= std::min(k, Ord1 - 1); j++)
      res[k] = res[k] + lhs[j] * rhs[k - j];
  return res;
}

///
/**
   Extract the 'comp'-th component.
   Require CoefType to support the [] syntax.
 */
template <int Order, class CoefType>
inline auto getComp(const Polynomial<Order, CoefType> &poly, int comp) {
  using T = std::decay_t<decltype(poly[0][0])>;  // the additional [0] is to get
                                                 // a component
  Polynomial<Order, T> result;
  for (int d = 0; d < Order; ++d) result[d] = poly[d][comp];
  return result;
}

//==================================================
// roots-finding functions for real polynomails

// template <int Order>
// using RealPolynomial = Polynomial<Order, Real>;

template <class CoefType, int Order>
struct _PolynomialRoots {
  template <class Inserter>
  static Inserter roots(const Polynomial<Order, CoefType> &, Inserter rts,
                        Real) {
    return rts;
  }
};

template <class CoefType>
struct _PolynomialRoots<CoefType, 2> {
  template <class Inserter>
  static Inserter roots(const Polynomial<2, CoefType> &poly, Inserter rts,
                        Real tol) {
    *rts++ = poly[0] / poly[1];
    return rts;
  }
};

template <class CoefType>
struct _PolynomialRoots<CoefType, 3> {
  template <class Inserter>
  static Inserter roots(const Polynomial<3, CoefType> &poly, Inserter rts,
                        Real tol) {
    const Real &a = poly[2];
    const Real &b = poly[1];
    const Real &c = poly[0];
    Real delta = b * b - 4 * a * c;
    if (delta < -tol) return rts;
    if (delta < tol) {
      *rts++ = -b / (2 * a);
      return rts;
    }
    Real dsqrt = sqrt(delta);
    *rts++ = (-b - dsqrt) / (2 * a);
    *rts++ = (-b + dsqrt) / (2 * a);
    return rts;
  }
};

template <class CoefType, class Inserter, int Order>
inline Inserter roots(const Polynomial<Order, CoefType> &poly, Inserter rts,
                      Real tol) {
  return _PolynomialRoots<CoefType, Order>::roots(poly, rts, tol);
};

template <class CoefType, int Order>
inline Real root(const Polynomial<Order, CoefType> &poly, Real x0, Real tol,
                 int maxIter) {
  auto q = poly.der();
  return fzero(poly, q, x0, maxIter, tol);
}

template <class CoefType, class Inserter, int Order>
inline Inserter extrema(const Polynomial<Order, CoefType> &poly, Inserter rts,
                        Real tol) {
  if (Order >= 3) return roots(poly.der(), rts, tol);
  return rts;
}

template <class T_Func, class T_Der>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
inline Real
fzero(const T_Func &hf, const T_Der &hdf, Real hx0, int maxIter, Real tol) {
  using localReal = Real;
  // Polynomial<decltype(hf)::_Order, localReal> f = hf;
  // Polynomial<decltype(hdf)::_Order, localReal> df = hdf;
  const auto &f = hf;
  const auto &df = hdf;
  localReal x0 = hx0;
  localReal fx = f(x0);
  localReal localFx = fx;
  localReal localX0 = x0;
  localReal dx = 0;
  int whileCount = 0;
  while (maxIter-- > 0 && (std::abs(fx) > tol)) {
    localReal dfx = df(x0);
    dx = fx / dfx;
    while (std::abs(localFx) >= std::abs(fx)) {
      localX0 = x0 - dx;
      localFx = f(localX0);
      if (whileCount++ > newtonSubCount())
        throw std::runtime_error(
            std::format("whileCount: Newton iteration may not converge, f(x) = "
                        "{}, dx = {} \n",
                        fx, dx));
      dx /= 2;
    }
    whileCount = 0;
    x0 = localX0;
    fx = localFx;
  }
  if (maxIter < 0 && (std::abs(fx) > tol))
    throw std::runtime_error(std::format(
        "maxIter: Newton iteration may not converge, f(x) = {}, dx = {} \n", fx,
        dx));
  return x0;
}

template <>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
inline Real
fzero<Polynomial<4, long double>, Polynomial<3, long double>>(
    const Polynomial<4, long double> &hf, const Polynomial<3, long double> &hdf,
    Real hx0, int maxIter, Real tol) {
  using localReal = long double;
  Polynomial<4, localReal> f(hf);
  auto df = f.der();
  localReal x0 = hx0;
  localReal fx = f(x0);
  localReal localFx = fx;
  localReal localX0 = x0;
  localReal dx = 0;
  localReal div = 1;
  int whileCount = 0;
  while (maxIter-- > 0 && (std::abs(fx) > tol)) {
    localReal dfx = df(x0);
    dx = fx / dfx;
    while (std::abs(localFx) >= std::abs(fx)) {
      localX0 = x0 - dx / div;
      localFx = f(localX0);
      div *= 2;
      if (whileCount++ > newtonSubCount()) {
        fzero<Polynomial<4, long double>, Polynomial<3, long double>>(
            hf, hdf, hx0, maxIter, tol);
        std::cout << std::format(
            "long double version, whileCount: Newton iteration may not "
            "converge, f(x) = "
            "{}, dx = {} \n",
            fx, dx / div);
        break;
      }
    }
    whileCount = 0;
    div /= 2;
    x0 = localX0;
    fx = localFx;
  }
  if (maxIter < 0 && (std::abs(fx) > tol))
    std::cout << std::format(
        "long double version, maxIter: Newton iteration may not converge, f(x) "
        "= {}, dx = {} \n",
        fx, dx / div);
  return x0;
}

// template <> struct _PolynomialRoots<Real, 2>;
// template <> struct _PolynomialRoots<long double, 2>;
// template <> struct _PolynomialRoots<Real, 3>;
// template <> struct _PolynomialRoots<long double, 3>;

#endif  // POLYNOMIAL_H
