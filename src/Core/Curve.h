#ifndef CURVE_H
#define CURVE_H

#include <vector>

#include "Core/Interval.h"
#include "Core/VecCompare.h"
#include "Polynomial.h"

template <int Dim, int Order>
class Curve {
 public:
  // fitCurve's boundary type
  enum BCType { notAKnot = 0, periodic, complete, second, nature, nBC_type };
  using rVec = Vec<Real, Dim>;
  using T_Polynomial = Polynomial<Order, rVec>;
  template <class T>
  using vector = std::vector<T>;

 protected:
  std::vector<Real> knots;
  std::vector<Polynomial<Order, rVec>> polys;

 public:
  // constructors
  Curve() = default;

  Curve(const Polynomial<Order, rVec>& pn, Real len) {
    polys.push_back(pn);
    knots.push_back(0);
    knots.push_back(len);
  }

  template<int Order2>
  explicit Curve(const Curve<Dim, Order2>&& rhs) : knots((rhs.getKnots())) {
    static_assert(Order2 < Order, "Order2 must be less than Order");
    for (const auto& p : rhs.getPolys()) {
      T_Polynomial poly(0);
      for (int i = 0; i < Order2; ++i) {
        poly[i] = p[i];
      }
      polys.push_back(poly);
    }
  }

  Curve(std::vector<Real>&& knots_,
        std::vector<Polynomial<Order, rVec>>&& polys_)
      : knots(std::move(knots_)), polys(std::move(polys_)) {}

  Curve(const Curve&) = default;
  auto operator=(const Curve&) -> Curve& = default;
  Curve(Curve&&) = default;
  auto operator=(Curve&&) -> Curve& = default;

  ~Curve() = default;

  // accessors
 public:
  // helpers
  int locatePiece(Real t) const;

OPTNONE_FUNC
  rVec
  operator()(Real t) const {
    int p = locatePiece(t);
    return (polys[p])(t - knots[p]);
  }

  Interval<1> domain() const {
    if (empty()) return Interval<1>(0.0, 0.0);
    return Interval<1>(knots.front(), knots.back());
  }

  bool empty() const { return knots.size() <= 1; }

  rVec startpoint() const { return polys.front()[0]; }

  rVec endpoint() const {
    auto np = polys.size();
    Real t = knots[np] - knots[np - 1];
    return polys[np - 1](t);
  }

  rVec midpoint() const {
    return (*this)(0.5 * (knots.front() + knots.back()));
  }

OPTNONE_FUNC
  bool
  isClosed(Real tol) const {
    auto poly = polys.back() - polys.front()[0];
    Real t = knots.back() - *(knots.end() - 2);
    auto res = poly(t);
    return VecCompare<Real, Dim>(tol).compare(res, rVec(0.0)) == 0;
    // VecCompare<Real, Dim> vCmp(tol);
    // return vCmp.compare(startpoint(), endpoint()) == 0;
  }

  const std::vector<Real>& getKnots() const { return knots; }
  auto getKnotPoints() const;

  const std::vector<T_Polynomial>& getPolys() const { return polys; }

  // modifiers
 public:
  void concat(const T_Polynomial& p, Real plen);

  void concat(const Curve<Dim, Order>& pp);

  Curve<Dim, Order> reverse() const;

  Curve<Dim, Order> offset(const rVec& ofs) const {
    Curve<Dim, Order> res;
    res.knots = this->knots;
    for (const auto& p : this->polys) res.polys.push_back(p + ofs);
    return res;
  }

  Curve<Dim, Order> rotate(Real theta) const {
    const Real c = cos(theta);
    const Real s = sin(theta);
    Curve<Dim, Order> res;
    res.knots = this->knots;
    for (const auto& p : this->polys) {
      T_Polynomial q;
      for (int k = 0; k < Order; ++k) {
        q[k][0] = p[k][0] * c + p[k][1] * (-s);
        q[k][1] = p[k][0] * s + p[k][1] * c;
      }
      res.polys.push_back(q);
    }
    return res;
  }

  // advanved operations
 public:
  // extact
  Curve<Dim, Order> extract(Real lo, Real hi, Real tol,
                            bool exact = false) const;

  void split(const vector<Real>& brks, vector<Curve<Dim, Order>>& out, Real tol,
             bool exact = false) const;
  ///
  /**
     Return a copy of this curve,
     so that each piece is monotonic in both x and y direction.
   */
  std::pair<Curve<Dim, Order>, vector<long double>> makeMonotonic(
      Real tol) const;
  ///
  /**
     Report the number of proper intersections with the line x_d=c.
     Require the curve to be piecewise monotonic, see makeMonotonic().
   */
  int countProperInts(Real c, int d, Real tol) const;

  // normal Vector
  std::vector<rVec> normalVector(Real brk, Real tol) const;

  // compare around a point p(boundary point)
  int compare(const Curve& rhs, const rVec& p, int ix, int iy, Real tol) const;
  // type = 0, startpoint()
  // type = 1, endpoint()
  auto getComparablePoint(Real tol, int type = 0) const -> rVec;
  bool equal(const Curve& rhs, Real tol) const;

  // curvature at t
  static Real curvature(const Polynomial<Order, Real>& xPoly,
                        const Polynomial<Order, Real>& yPoly, Real t);

 protected:
  static Real paraCalculator(const auto& poly, int ix, Real x, Real t0,
                             Real tol) {
    auto xPoly = getComp(poly, ix);
    return root(xPoly - x, t0, tol, newtonMaxIter());
  };

  // friend helpers
 public:
  template <int Ord>
  friend Curve<2, Ord> createRect(const Vec<Real, 2>& lo,
                                  const Vec<Real, 2>& hi);
  template <int Ord>
  friend Curve<2, Ord> createLineSegment(const Vec<Real, 2>& p0,
                                         const Vec<Real, 2>& p1);
  template <int Dm, int Ord>
  friend Curve<Dm, Ord - 1> der(const Curve<Dm, Ord>& c);
  // template <int Ord>
  // friend Curve<2,Ord> fitCurve(const std::vector<Vec<Real,2>> &knots,
  //                              bool periodic);
  template <int Ord>
  friend Curve<2, Ord> fitCurve(const std::vector<Vec<Real, 2>>& knots,
                                typename Curve<2, Ord>::BCType type,
                                const Vec<Real, 2>& start,
                                const Vec<Real, 2>& end);
  // read/write operations
 public:
  static Curve<Dim, Order> load(std::istream& is);
  void dump(std::ostream& os) const;
};

//================================================

template <int Dim, int Order>
Curve<Dim, Order - 1> der(const Curve<Dim, Order>& c);

template <int Dim, int Order>
Interval<Dim> boundingBox(const Curve<Dim, Order>& c);

template <int Dim, int Order>
Interval<Dim> boundingBox(const std::vector<Curve<Dim, Order>>& vc);

template <int Order>
Real area(const Curve<2, Order>& gon, Real tol = distTol());

template <int Order>
Real arclength(const Curve<2, Order>& c);

template <int Ord>
Curve<2, Ord> createRect(const Vec<Real, 2>& lo, const Vec<Real, 2>& hi);

template <int Ord>
Curve<2, Ord> createLineSegment(const Vec<Real, 2>& p0, const Vec<Real, 2>& p1);

template <int Order>
Curve<2, Order> fitCurve(
    const std::vector<Vec<Real, 2>>& knots,
    typename Curve<2, Order>::BCType type = Curve<2, Order>::notAKnot,
    const Vec<Real, 2>& start = Vec<Real, 2>(),
    const Vec<Real, 2>& end = Vec<Real, 2>());

template <int Order>
void checkFitCurve(const Curve<2, Order>& crv,
                   typename Curve<2, Order>::BCType type,
                   const Vec<Real, 2>& start = Vec<Real, 2>(),
                   const Vec<Real, 2>& end = Vec<Real, 2>());

template <int Dim, int Order>
auto Curve<Dim, Order>::getKnotPoints() const {
  std::vector<rVec> ret;
  ret.reserve(knots.size());
  for (auto& poly : polys) {
    ret.push_back(poly[0]);
  }
  ret.push_back(endpoint());
  return ret;
}

#endif  // CURVE_H
