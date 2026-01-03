#include "Curve.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "Core/HighPrecisionAdd.h"
#include "Core/SolveTri.h"
#include "Core/VecCompare.h"
// #include "Wrapper_LAPACKE.h"

template <int Dim, int Order>
int Curve<Dim, Order>::locatePiece(Real t) const {
  int np = polys.size();
  int left = 0;
  int right = np - 1;
  int mid;
  for (;;) {
    if (left == right) return left;
    mid = (left + right) / 2;
    if (t < knots[mid + 1])
      right = mid;
    else
      left = mid + 1;
  }
}

template <int Dim, int Order>
void Curve<Dim, Order>::concat(const T_Polynomial &p, Real plen) {
  if (knots.empty()) knots.push_back(0.0);
  polys.push_back(p);
  knots.push_back(knots.back() + plen);
}

template <int Dim, int Order>
void Curve<Dim, Order>::concat(const Curve<Dim, Order> &pp) {
  if (pp.empty()) return;
  if (knots.empty()) knots.push_back(0.0);
  const auto &newknots = pp.getKnots();
  const auto &newpolys = pp.getPolys();
  Real delta = knots.back() - newknots.front();
  for (std::size_t i = 0; i < newpolys.size(); ++i) {
    polys.push_back(newpolys[i]);
    knots.push_back(newknots[i + 1] + delta);
  }
}

template <int Dim, int Order>
Curve<Dim, Order> Curve<Dim, Order>::reverse() const {
  Curve<Dim, Order> res;
  res.knots.resize(knots.size());
  res.polys.resize(polys.size());

  Real lastknot = knots.back();
  int j = polys.size() - 1;
  res.knots[0] = 0;
  for (std::size_t i = 0; i < polys.size(); i++, j--) {
    // flip the knots
    res.knots[i + 1] = lastknot - knots[j];
    // flip the parametrization
    res.polys[i] =
        polys[j].translate(knots[j + 1] - knots[j], true);  // reflect = true
  }
  return res;
}

/*
template <int Dim, int Order>
Curve<Dim, Order>
Curve<Dim, Order>::extract(Real lo, Real hi, Real tol) const
{
  if (hi <= lo + tol)
    return Curve<Dim, Order>();

  int ihead = locatePiece(lo);
  int itail = locatePiece(hi);
  // avoid pieces with tiny length
  if (lo >= knots[ihead + 1] - tol)
    ihead++;
  if (hi <= knots[itail] + tol)
    itail--;
  if (ihead > itail)
  {
    std::cout << knots[ihead] << std::endl;
    std::cout << knots[itail] << std::endl;
    assert(ihead <= itail);
  }
  //assert(ihead <= itail);

  Curve<Dim, Order> res;
  res.knots.resize(itail - ihead + 2);
  res.polys.resize(itail - ihead + 1);
  res.knots[0] = lo;
  res.polys[0] = polys[ihead].translate(lo - knots[ihead]);

  // range [ihead+1, itail]
  std::copy(&knots[ihead + 1], &knots[itail + 1], &(res.knots[1]));
  std::copy(&polys[ihead + 1], &polys[itail + 1], &(res.polys[1]));
  // mark the tail
  res.knots.back() = hi;

  return res;
}
*/

template <int Dim, int Order>
OPTNONE_FUNC
Curve<Dim, Order>
Curve<Dim, Order>::extract(Real lo, Real hi, Real tol, bool exact) const {
  if (hi <= lo) return Curve<Dim, Order>();
  if (!exact && hi <= lo + tol) return Curve<Dim, Order>();

  int ihead = locatePiece(lo);
  int itail = locatePiece(hi);

  // if ((lo >= knots[ihead + 1] - tol) && (hi <= knots[itail] + tol) &&
  //     (ihead + 1 == itail)) {
  //   Curve<Dim, Order> res;
  //   res.knots.resize(2);
  //   res.polys.resize(1);
  //   res.knots[0] = lo;
  //   res.polys[0] = polys[ihead + 1].translate(lo - knots[ihead + 1]);
  //   res.knots[1] = hi;
  //   return res;
  // }

  // avoid pieces with tiny length
  if (!exact && lo >= knots[ihead + 1] - tol) ihead++;
  if (!exact && hi <= knots[itail] + tol) itail--;
  if (ihead > itail) {
    throw std::runtime_error("Curve::extract: invalid range.");
  }
  // assert(ihead <= itail);

  Curve<Dim, Order> res;
  res.knots.resize(itail - ihead + 2);
  res.polys.resize(itail - ihead + 1);
  res.knots[0] = lo;
  res.polys[0] = polys[ihead].translate(lo - knots[ihead]);

  // range [ihead+1, itail]
  std::copy(&knots[ihead + 1], &knots[itail + 1], &(res.knots[1]));
  std::copy(&polys[ihead + 1], &polys[itail + 1], &(res.polys[1]));
  // mark the tail
  res.knots.back() = hi;

  return res;
}

template <int Dim, int Order>
OPTNONE_FUNC
inline void
Curve<Dim, Order>::split(const vector<Real> &brks,
                         vector<Curve<Dim, Order>> &out, Real tol,
                         bool exact) const {
  bool closed = isClosed(tol);
  if (exact) closed = isClosed(distTol());
  if (brks.empty()) {
    out.push_back(*this);
    return;
  }
  auto head = extract(knots.front(), brks.front(), tol, exact);
  if (!closed && !head.empty()) out.push_back(std::move(head));
  std::size_t i = 0;
  for (; i < brks.size() - 1; ++i) {
    if (exact || brks[i + 1] > brks[i] + tol)
      out.push_back(extract(brks[i], brks[i + 1], tol, exact));
  }
  auto tail = extract(brks[i], knots.back(), tol, exact);
  if (!closed) {
    if (!tail.empty()) out.push_back(std::move(tail));
  } else {
    tail.concat(head);
    if (!tail.empty()) out.push_back(std::move(tail));
  }
}

template <int Dim, int Order>
OPTNONE_FUNC
std::pair<Curve<Dim, Order>, std::vector<long double>>
Curve<Dim, Order>::makeMonotonic(Real tol) const {
  using localReal = long double;
  Curve<Dim, Order> res;
  std::vector<localReal> newKnotsDist; // for high precision
  int np = knots.size() - 1;
  for (int i = 0; i < np; i++) {
    // auto p0 = polys[i](0);
    // auto p1 = polys[i](knots[i + 1] - knots[i]);
    Polynomial<Order, Vec<localReal, Dim>> polyi(polys[i]);
    std::vector<localReal> ex;
    // find out all the monotonic pieces by first locating the extrema
    for (int d = 0; d < Dim; d++) {
      auto rp = getComp(polyi, d);
      extrema<localReal>(rp, std::back_inserter(ex), tol);
      // std::vector<localReal> localEx;
      // extrema<localReal>(rp, std::back_inserter(localEx), tol);
      // localReal bound[2] = {rp(0), rp(knots[i + 1] - knots[i])};
      // if (bound[0] > bound[1]) std::swap(bound[0], bound[1]);
      // for (auto& t : localEx) {
      //   auto v = rp(t);
      //   if (v < bound[0] - tol || v > bound[1] + tol) ex.push_back(t);
      // }
    }
    // filter out the extrema out of domain
    // note that polys[i] is expressed in the variable (t-knots[i])
    Interval<1> validDomain(0, knots[i + 1] - knots[i]);
    auto tbr1 =
        std::remove_if(ex.begin(), ex.end(), [&](const Real &a) -> bool {
          return !validDomain.contain(a, (Real)0);
        });
    // zero-tolerance is fine ^
    std::sort(ex.begin(), tbr1);
    ex.erase(std::unique(
                 ex.begin(), tbr1,
                 [&tol](const Real &a, const Real &b) { return b - a < tol; }),
             ex.end());

    if (!ex.empty()) {
      // mark the head & tail while avoiding tiny pieces
      if (ex.front() > tol)
        ex.insert(ex.begin(), (Real)0);
      else
        ex.front() = 0;
      if (ex.back() < validDomain.hi()[0] - tol)
        ex.push_back(validDomain.hi()[0]);
      else
        ex.back() = validDomain.hi()[0];
      for (std::size_t j = 0; j < ex.size() - 1; j++) {
        T_Polynomial polyj(polyi.translate(ex[j]));
        res.concat(polyj, ex[j + 1] - ex[j]);
        newKnotsDist.push_back(ex[j + 1] - ex[j]);
      }
    } else {
      // just copy if it is already mono
      res.concat(polys[i], knots[i + 1] - knots[i]);
      newKnotsDist.push_back(knots[i + 1] - knots[i]);
    }
  }
  return {res, newKnotsDist};
}

template <int Dim, int Order>
int Curve<Dim, Order>::countProperInts(Real c, int d, Real tol) const {
  int res = 0;
  auto process = [&c, &d, &tol](const rVec &p0, const rVec &p1) -> int {
    Real lo = std::min(p0[d], p1[d]);
    Real hi = std::max(p0[d], p1[d]);
    if (hi - lo <= tol)  // parallel to axis
      return 0;
    if (hi < c - tol || lo > c + tol)  // disjoint case
      return 0;
    if (hi <= c + tol || lo >= c - tol) {  // improper case
      auto mid = (p0 + p1) * 0.5;
      return (mid[d] > c) ? (1) : (0);
    }
    return 1;
  };
  for (std::size_t i = 0; i < polys.size(); ++i)
    res += process(polys[i][0], polys[i](knots[i + 1] - knots[i]));
  return res;
}

//============================================================

template <int Dim, int Order>
void Curve<Dim, Order>::dump(std::ostream &os) const {
  int header[] = {Dim, Order};
  os.write((char *)header, sizeof(header));

  int np = polys.size();
  os.write((char *)&np, sizeof(np));

  os.write((char *)knots.data(), sizeof(Real) * (np + 1));
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < Dim; j++) {
      auto rp = getComp(polys[i], j);
      os.write((char *)(&rp[0]), sizeof(Real) * Order);
    }
  }
}

template <int Dim, int Order>
Curve<Dim, Order> Curve<Dim, Order>::load(std::istream &is) {
  Curve<Dim, Order> res;
  std::vector<Real> &knots = res.knots;
  std::vector<Polynomial<Order, rVec>> &polys = res.polys;

  int tmp[2];
  int np;
  is.read((char *)tmp, sizeof(int) * 2);
  assert(tmp[0] == Dim && tmp[1] == Order);
  is.read((char *)&np, sizeof(int));
  knots.resize(np + 1);
  polys.resize(np);
  is.read((char *)knots.data(), sizeof(Real) * (np + 1));
  for (int i = 0; i < np; i++) {
    Polynomial<Order, rVec> p;
    Real buf[Dim * Order];
    Real *pbuf = buf;
    is.read((char *)buf, sizeof(Real) * Dim * Order);
    for (int j = 0; j < Dim; j++)
      for (int k = 0; k < Order; k++) p[k][j] = *pbuf++;
    polys[i] = p;
  }
  return res;
}

//============================================================

template <int Order>
Curve<2, Order> createLineSegment(const Vec<Real, 2> &p0,
                                  const Vec<Real, 2> &p1) {
  Curve<2, Order> res;
  Real l = norm(p1 - p0, 2);
  res.knots = {0.0, l};
  res.polys.resize(1);
  res.polys.front()[0] = p0;
  res.polys.front()[1] = (p1 - p0) / l;
  return res;
};

template <int Order>
Curve<2, Order> createRect(const Vec<Real, 2> &lo, const Vec<Real, 2> &hi) {
  Curve<2, Order> res;
  auto &knots = res.knots;
  auto &polys = res.polys;
  res.knots.resize(5);
  res.polys.resize(4);
  knots[0] = 0.0;
  knots[1] = hi[0] - lo[0];
  knots[2] = knots[1] + hi[1] - lo[1];
  knots[3] = knots[2] + knots[1];
  knots[4] = 2 * knots[2];
  polys[0][0] = lo;
  polys[0][1] = {1.0, 0.0};
  polys[1][0] = {hi[0], lo[1]};
  polys[1][1] = {0.0, 1.0};
  polys[2][0] = hi;
  polys[2][1] = {-1.0, 0.0};
  polys[3][0] = {lo[0], hi[1]};
  polys[3][1] = {0.0, -1.0};
  return res;
}

template <int Order>
OPTNONE_FUNC
Real area(const Curve<2, Order> &gon, Real tol) {
  Real a = 0.0;
  if (gon.empty()) return a;
  if (!gon.isClosed(tol)) {
    gon.isClosed(tol);
    throw std::runtime_error("unclosed Curve calculate area. may enlarge Tol value.");
  }
  const auto &knots = gon.getKnots();
  auto pts = gon.getKnotPoints();
  Real xl = pts.front()[0];
  Real xh = pts.front()[0];
  Real yl = pts.front()[1];
  Real yh = pts.front()[1];
  for (auto &p : pts) {
    xl = std::min(xl, p[0]);
    xh = std::max(xh, p[0]);
    yl = std::min(yl, p[1]);
    yh = std::max(yh, p[1]);
  }
  bool xdy = true;
  if (xh - xl > yh - yl) xdy = false;
  int i = 0;
  std::vector<Real> vals;
  vals.reserve(gon.getPolys().size());
  // apply the Green's formula
  for (const auto &p : gon.getPolys()) {
    auto x = getComp(p, 0);
    auto y = getComp(p, 1);
    Real localVal;
    if (xdy) {
      auto dy = y.der(); 
      localVal = (x * dy).prim()(knots[i + 1] - knots[i]);
    } else {
      auto dx = x.der();
      localVal = (-y * dx).prim()(knots[i + 1] - knots[i]);
    }
    // a += localVal;
    vals.push_back(localVal);
    ++i;
  }
  a = highPrecisionAdd(vals);
  return a;
}

template <int Order>
Real arclength(const Curve<2, Order> &c) {
  Real l = 0.0;
  if (c.empty()) return l;
  const auto &knots = c.getKnots();
  const auto &polys = c.getPolys();
  for (std::size_t k = 0; k < polys.size(); ++k) {
    auto dt = knots[k + 1] - knots[k];
    auto dxdt = getComp(polys[k], 0).der();
    auto dydt = getComp(polys[k], 1).der();
    auto ndsdt = [&](Real t) {
      Vec<Real, 2> ds{dxdt(t), dydt(t)};
      return norm(ds, 2);
    };
    l += aquad(ndsdt, 0, dt);
  }
  return l;
}

template <int Dim, int Order>
Interval<Dim> boundingBox(const Curve<Dim, Order> &c) {
  std::vector<Curve<Dim, Order>> vc{c};
  return boundingBox(vc);
}

template <int Dim, int Order>
Interval<Dim> boundingBox(const std::vector<Curve<Dim, Order>> &vc) {
  Vec<Real, Dim> lower(std::numeric_limits<Real>::max());
  Vec<Real, Dim> upper = -lower;
  for (const auto &c : vc) {
    const auto &polys = c.getPolys();
    for (std::size_t i = 0; i < polys.size(); ++i) {
      lower = min(lower, polys[i][0]);
      upper = max(upper, polys[i][0]);
    }
    lower = min(lower, c.endpoint());
    upper = max(upper, c.endpoint());
  }
  return Interval<Dim>(lower, upper);
}

template <int Dim, int Order>
Curve<Dim, Order - 1> der(const Curve<Dim, Order> &c) {
  Curve<Dim, Order - 1> res;
  res.knots = c.getKnots();
  for (const auto &pn : c.getPolys()) res.polys.push_back(pn.der());
  return res;
}

template <>
Curve<2, 2> fitCurve(const std::vector<Vec<Real, 2>> &knots,
                     typename Curve<2, 2>::BCType, const Vec<Real, 2> &,
                     const Vec<Real, 2> &) {
  const int Order = 2;
  auto numKnots = knots.size();
  assert(numKnots >= 2);
  Curve<2, Order> res;
  for (std::size_t i = 0; i < numKnots - 1; ++i) {
    Polynomial<2, Vec<Real, 2>> p;
    Real l = norm(knots[i + 1] - knots[i], 2);
    p[0] = knots[i];
    p[1] = (knots[i + 1] - knots[i]) / l;
    res.concat(p, l);
  }
  return res;
}

template <>
Curve<2, 4> fitCurve(const std::vector<Vec<Real, 2>> &vertices,
                     typename Curve<2, 4>::BCType type,
                     const Vec<Real, 2> &start, const Vec<Real, 2> &end) {
  // using Point = Vec<Real, 2>;
  const int Order = 4;
  const int numPiece = vertices.size() - 1;
  int k;
  // calculate the accumulated chordal length
  std::vector<Real> t(numPiece + 1);
  t[0] = 0.0;
  for (k = 1; k <= numPiece; ++k)
    t[k] = t[k - 1] + norm(vertices[k] - vertices[k - 1], 2);
  if (type == Curve<2, 4>::periodic) {
    // prepare the coefficient matrix
    std::vector<Real> a(numPiece), b(numPiece, 2.0), c(numPiece);
    a[0] = t[1] / (t[1] + t[numPiece] - t[numPiece - 1]);
    c[0] = (t[numPiece] - t[numPiece - 1]) /
           (t[1] + t[numPiece] - t[numPiece - 1]);
    for (k = 1; k < numPiece; ++k) {
      a[k] = (t[k + 1] - t[k]) / (t[k + 1] - t[k - 1]);
      c[k] = (t[k] - t[k - 1]) / (t[k + 1] - t[k - 1]);
    }
    // prepare the RHS
    std::vector<Real> rhsx(numPiece), rhsy(numPiece);
    const auto &verts = vertices;
    rhsx[0] = 3 *
                  ((t[numPiece] - t[numPiece - 1]) /
                   (t[1] + t[numPiece] - t[numPiece - 1])) *
                  ((verts[1][0] - verts[0][0]) / t[1]) +
              3 * (t[1] / (t[1] + t[numPiece] - t[numPiece - 1])) *
                  ((verts[0][0] - verts[numPiece - 1][0]) /
                   (t[numPiece] - t[numPiece - 1]));
    rhsy[0] = 3 *
                  ((t[numPiece] - t[numPiece - 1]) /
                   (t[1] + t[numPiece] - t[numPiece - 1])) *
                  ((verts[1][1] - verts[0][1]) / t[1]) +
              3 * (t[1] / (t[1] + t[numPiece] - t[numPiece - 1])) *
                  ((verts[0][1] - verts[numPiece - 1][1]) /
                   (t[numPiece] - t[numPiece - 1]));
    for (int k = 1; k < numPiece; ++k) {
      rhsx[k] = 3 * ((t[k] - t[k - 1]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k + 1][0] - verts[k][0]) / (t[k + 1] - t[k])) +
                3 * ((t[k + 1] - t[k]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k][0] - verts[k - 1][0]) / (t[k] - t[k - 1]));
      rhsy[k] = 3 * ((t[k] - t[k - 1]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k + 1][1] - verts[k][1]) / (t[k + 1] - t[k])) +
                3 * ((t[k + 1] - t[k]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k][1] - verts[k - 1][1]) / (t[k] - t[k - 1]));
    }
    // solve the linear system :use LUD decomposition to solve the
    // periodic-tridiagonal system
    std::vector<Real> resx = solvePeriodicTri(a, b, c, rhsx);
    std::vector<Real> resy = solvePeriodicTri(a, b, c, rhsy);

    // assemble the spline
    Curve<2, Order> res;
    res.knots = std::vector<Real>(t.begin(), t.end());
    res.polys.resize(numPiece);
    for (k = 0; k < numPiece - 1; ++k) {
      Vec<Real, 2> K({(verts[k + 1][0] - verts[k][0]) / (t[k + 1] - t[k]),
                      (verts[k + 1][1] - verts[k][1]) / (t[k + 1] - t[k])});
      Vec<Real, 2> m1({resx[k], resy[k]});
      Vec<Real, 2> m2({resx[k + 1], resy[k + 1]});
      res.polys[k][3] =
          (m1 + m2 - (K * 2)) / ((t[k + 1] - t[k]) * (t[k + 1] - t[k]));
      res.polys[k][2] = (K * 3 - (m1 * 2) - m2) / (t[k + 1] - t[k]);
      res.polys[k][1] = m1;
      res.polys[k][0] = Vec<Real, 2>{verts[k][0], verts[k][1]};
    }
    Vec<Real, 2> K({(verts[numPiece][0] - verts[numPiece - 1][0]) /
                        (t[numPiece] - t[numPiece - 1]),
                    (verts[numPiece][1] - verts[numPiece - 1][1]) /
                        (t[numPiece] - t[numPiece - 1])});
    Vec<Real, 2> m1({resx[numPiece - 1], resy[numPiece - 1]});
    Vec<Real, 2> m2({resx[0], resy[0]});
    res.polys[numPiece - 1][3] =
        (m1 + m2 - (K * 2)) /
        ((t[numPiece] - t[numPiece - 1]) * (t[numPiece] - t[numPiece - 1]));
    res.polys[numPiece - 1][2] =
        (K * 3 - (m1 * 2) - m2) / (t[numPiece] - t[numPiece - 1]);
    res.polys[numPiece - 1][1] = m1;
    res.polys[numPiece - 1][0] =
        Vec<Real, 2>{verts[numPiece - 1][0], verts[numPiece - 1][1]};
    return res;
  } else if (type == Curve<2, 4>::notAKnot || type == Curve<2, 4>::complete ||
             type == Curve<2, 4>::second || type == Curve<2, 4>::nature) {
    // prepare the coefficient matrix
    std::vector<Real> a(numPiece), b(numPiece + 1, 2.0), c(numPiece);
    Real d = 0, e = 0, tmp;
    if (type == Curve<2, 4>::notAKnot) {
      if (numPiece <= 2)
        throw std::runtime_error("notAknot curve must have at least 3 pieces");
      tmp = (t[1] - t[0]) * (t[1] - t[0]) / ((t[2] - t[1]) * (t[2] - t[1]));
      b[0] = 1.0;
      c[0] = 1.0 - tmp;
      d = -tmp;
      tmp = (t[numPiece] - t[numPiece - 1]) * (t[numPiece] - t[numPiece - 1]) /
            ((t[numPiece - 1] - t[numPiece - 2]) *
             (t[numPiece - 1] - t[numPiece - 2]));
      b[numPiece] = 1.0;
      a[numPiece - 1] = 1.0 - tmp;
      e = -tmp;
    } else if (type == Curve<2, 4>::complete) {
      b[0] = 1.0;
      c[0] = 0.0;
      b[numPiece] = 1.0;
      a[numPiece - 1] = 0.0;
    } else if (type == Curve<2, 4>::second || type == Curve<2, 4>::nature) {
      c[0] = 1.0;
      a[numPiece - 1] = 1.0;
    }
    for (k = 1; k < numPiece; ++k) {
      a[k - 1] = (t[k + 1] - t[k]) / (t[k + 1] - t[k - 1]);
      c[k] = (t[k] - t[k - 1]) / (t[k + 1] - t[k - 1]);
    }
    // prepare the RHS
    std::vector<Real> rhsx(numPiece + 1), rhsy(numPiece + 1);
    const auto &verts = vertices;
    if (type == Curve<2, 4>::notAKnot) {
      tmp = (t[1] - t[0]) * (t[1] - t[0]) / ((t[2] - t[1]) * (t[2] - t[1]));
      rhsx[0] = -2 * tmp * ((verts[2][0] - verts[1][0]) / (t[2] - t[1])) +
                2 * ((verts[1][0] - verts[0][0]) / (t[1] - t[0]));
      rhsy[0] = -2 * tmp * ((verts[2][1] - verts[1][1]) / (t[2] - t[1])) +
                2 * ((verts[1][1] - verts[0][1]) / (t[1] - t[0]));
      tmp = (t[numPiece] - t[numPiece - 1]) * (t[numPiece] - t[numPiece - 1]) /
            ((t[numPiece - 1] - t[numPiece - 2]) *
             (t[numPiece - 1] - t[numPiece - 2]));
      rhsx[numPiece] = -2 * tmp *
                           ((verts[numPiece - 1][0] - verts[numPiece - 2][0]) /
                            (t[numPiece - 1] - t[numPiece - 2])) +
                       2 * ((verts[numPiece][0] - verts[numPiece - 1][0]) /
                            (t[numPiece] - t[numPiece - 1]));
      rhsy[numPiece] = -2 * tmp *
                           ((verts[numPiece - 1][1] - verts[numPiece - 2][1]) /
                            (t[numPiece - 1] - t[numPiece - 2])) +
                       2 * ((verts[numPiece][1] - verts[numPiece - 1][1]) /
                            (t[numPiece] - t[numPiece - 1]));
    } else if (type == Curve<2, 4>::complete) {
      rhsx[0] = start[0];
      rhsy[0] = start[1];
      rhsx[numPiece] = end[0];
      rhsy[numPiece] = end[1];
    } else if (type == Curve<2, 4>::nature) {
      rhsx[0] = 3 * ((verts[1][0] - verts[0][0]) / (t[1] - t[0]));
      rhsy[0] = 3 * ((verts[1][1] - verts[0][1]) / (t[1] - t[0]));
      rhsx[numPiece] = 3 * ((verts[numPiece][0] - verts[numPiece - 1][0]) /
                            (t[numPiece] - t[numPiece - 1]));
      rhsy[numPiece] = 3 * ((verts[numPiece][1] - verts[numPiece - 1][1]) /
                            (t[numPiece] - t[numPiece - 1]));
    } else if (type == Curve<2, 4>::second) {
      rhsx[0] = 3 * ((verts[1][0] - verts[0][0]) / (t[1] - t[0])) -
                0.5 * start[0] * (t[1] - t[0]);
      rhsy[0] = 3 * ((verts[1][1] - verts[0][1]) / (t[1] - t[0])) -
                0.5 * start[1] * (t[1] - t[0]);
      rhsx[numPiece] = 3 * ((verts[numPiece][0] - verts[numPiece - 1][0]) /
                            (t[numPiece] - t[numPiece - 1])) -
                       0.5 * end[0] * (t[numPiece] - t[numPiece - 1]);
      rhsy[numPiece] = 3 * ((verts[numPiece][1] - verts[numPiece - 1][1]) /
                            (t[numPiece] - t[numPiece - 1])) -
                       0.5 * end[1] * (t[numPiece] - t[numPiece - 1]);
    }
    for (int k = 1; k < numPiece; ++k) {
      rhsx[k] = 3 * ((t[k] - t[k - 1]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k + 1][0] - verts[k][0]) / (t[k + 1] - t[k])) +
                3 * ((t[k + 1] - t[k]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k][0] - verts[k - 1][0]) / (t[k] - t[k - 1]));
      rhsy[k] = 3 * ((t[k] - t[k - 1]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k + 1][1] - verts[k][1]) / (t[k + 1] - t[k])) +
                3 * ((t[k + 1] - t[k]) / (t[k + 1] - t[k - 1])) *
                    ((verts[k][1] - verts[k - 1][1]) / (t[k] - t[k - 1]));
    }
    // solve the linear system
    std::vector<Real> resx;
    std::vector<Real> resy;
    if (type == Curve<2, 4>::notAKnot) {
      resx = solveTrisp(a, b, c, d, e, rhsx);
      resy = solveTrisp(a, b, c, d, e, rhsy);
    } else {
      resx = solveTri(a, b, c, rhsx);
      resy = solveTri(a, b, c, rhsy);
    }

    // assemble the spline
    Curve<2, Order> res;
    res.knots = std::vector<Real>(t.begin(), t.end());
    res.polys.resize(numPiece);
    for (k = 0; k < numPiece; ++k) {
      Vec<Real, 2> K({(verts[k + 1][0] - verts[k][0]) / (t[k + 1] - t[k]),
                      (verts[k + 1][1] - verts[k][1]) / (t[k + 1] - t[k])});
      Vec<Real, 2> m1({resx[k], resy[k]});
      Vec<Real, 2> m2({resx[k + 1], resy[k + 1]});
      res.polys[k][3] =
          (m1 + m2 - (K * 2)) / ((t[k + 1] - t[k]) * (t[k + 1] - t[k]));
      res.polys[k][2] = (K * 3 - (m1 * 2) - m2) / (t[k + 1] - t[k]);
      res.polys[k][1] = m1;
      res.polys[k][0] = Vec<Real, 2>{verts[k][0], verts[k][1]};
    }
    return res;
  } else
    exit(-1);
}

template <int Dim, int Order>
auto Curve<Dim, Order>::normalVector(Real brk, Real tol) const
    -> std::vector<Vec<Real, Dim>> {
  auto i = locatePiece(brk);
  auto &poly = polys[i];
  auto t = brk - knots[i];
  auto dx = getComp(poly, 0).der();
  auto dy = getComp(poly, 1).der();

  rVec ret = {dy(t), -dx(t)};
  return {ret / norm(ret)};
}

template <int Dim, int Order>
int Curve<Dim, Order>::compare(const Curve &rhs, const rVec &p, int ix, int iy,
                               Real tol) const {
  // extract the curve piece around the point p.
  VecCompare<Real, Dim> vcmp(tol);
  auto poly = polys[0];
  Real lo = knots[0];
  Real hi = knots[1];
  Real range = hi - lo;
  Real t0;
  if (polys.size() > 1) {
    if (vcmp(p, endpoint()) == 0) {
      poly = polys.back();
      lo = knots[knots.size() - 2];
      hi = knots.back();
    } else if (vcmp(p, startpoint()) != 0) {
      assert(false && "point p not curve boundary point.");
    }
  }
  if (vcmp(poly[0], p) == 0)
    t0 = 0;
  else
    t0 = range;

  if (rhs.getPolys().size() > 1)
    return -rhs.compare(Curve(poly, hi - lo), p, ix, iy, tol);

  auto rhsPoly = rhs.getPolys()[0];
  auto rhsLo = rhs.getKnots()[0];
  auto rhsHi = rhs.getKnots()[1];
  auto rhsRange = rhsHi - rhsLo;
  Real rhsT0;
  if (vcmp(rhsPoly[0], p) == 0)
    rhsT0 = 0;
  else
    rhsT0 = rhsRange;

  // find intersection with x = p[ix] - tol / 2.
  auto t = paraCalculator(poly, ix, p[ix] - tol / 2, t0, tol);
  auto rhsT = paraCalculator(rhsPoly, ix, p[ix] - tol / 2, rhsT0, tol);
  if ((t > 0 && t < range) && (rhsT > 0 && rhsT < rhsRange))
    return poly(t)[iy] < rhsPoly(rhsT)[iy] ? -1 : 1;
  // find intersection with x = p[ix] + tol / 2.
  t = paraCalculator(poly, ix, p[ix] + tol / 2, t0, tol);
  rhsT = paraCalculator(rhsPoly, ix, p[ix] + tol / 2, rhsT0, tol);
  if ((t > 0 && t < range) && (rhsT > 0 && rhsT < rhsRange))
    return poly(t)[iy] < rhsPoly(rhsT)[iy] ? -1 : 1;

  // uncomparable in (ix, iy) direction.
  return 0;
}

template <int Dim, int Order>
OPTNONE_FUNC
auto Curve<Dim, Order>::getComparablePoint(Real tol, int type) const
    -> rVec {
  if (type == 0) {
    int i = 0;
    // ignore tiny curve pieces
    while (i < polys.size() - 1 && tol > knots[i + 1] - knots[0]) ++i;
    Real lastRange = knots[i + 1] - knots[i];
    Real disturbance = std::min(tol, lastRange) / 2;
    Real t0 = disturbance;
    return normalize(polys[i](t0));
  }
  if (type == 1) {
    int i = polys.size() - 1;
    while (i > 0 && tol > knots[knots.size() - 1] - knots[i]) --i;
    Real lastRange = knots[i + 1] - knots[i];
    Real disturbance = std::min(tol, lastRange) / 2;
    Real t0 = lastRange - disturbance;
    return normalize(polys[i](t0));
  }
  throw std::runtime_error("type error");
  return rVec();
}

template <int Dim, int Order>
OPTNONE_FUNC
bool Curve<Dim, Order>::equal(const Curve &rhs, Real tol) const {
  // auto d1 = this->getComparablePoint(tol, 0);
  // auto d2 = rhs.getComparablePoint(tol, 0);
  // return (norm(d1 - d2) < tol);
  VecCompare<Real, Dim> vCmp(tol);
  if (vCmp.compare(midpoint(), rhs.midpoint()) != 0) return false;
  const auto &lhs = *this;
  auto &lhsKnots = this->getKnots();
  auto &rhsKnots = rhs.getKnots();
  const int numPoint = std::max(3, Order);  // at least contain a midpoint.
  Real lhsDt = (lhsKnots.back() - lhsKnots.front()) / (numPoint - 1);
  Real rhsDt = (rhsKnots.back() - rhsKnots.front()) / (numPoint - 1);
  for (int i = 0; i < numPoint; ++i) {
    auto pl = lhs(lhsKnots.front() + i * lhsDt);
    auto pr = rhs(rhsKnots.front() + i * rhsDt);
    if (vCmp.compare(pl, pr) != 0) return false;
  }
  return true;
}

template <int Dim, int Order>
auto Curve<Dim, Order>::curvature(const vector<Real> &params) const
    -> vector<Real> {
  vector<Real> res;
  if (params.empty() || empty()) return res;
  res.reserve(params.size());
  for (const auto &t : params) {
    const int piece = locatePiece(t);
    const auto &poly = polys[piece];
    const Real localT = t - knots[piece];
    res.push_back(
        Curve<Dim, Order>::curvature(getComp(poly, 0), getComp(poly, 1),
                                     localT));
  }
  return res;
}

namespace {
template <int Dim, int Order>
static bool insertByUpper(const Curve<Dim, Order> &curve,
                          const std::vector<Real> &upperBound,
                          std::vector<Real> &params,
                          std::vector<int> &fixed, Real tol) {
  if (params.size() < 2) return false;

  struct Insertion {
    std::size_t index;
    std::vector<Real> mids;
  };
  std::vector<Insertion> pending;
  pending.reserve(params.size());

  auto segDist = [&](std::size_t i, std::size_t j) {
    return norm(curve(params[j]) - curve(params[i]), 2);
  };

  for (std::size_t i = 0; i + 1 < params.size(); ++i) {
    const Real upperSeg = std::min(upperBound[i], upperBound[i + 1]);
    if (upperSeg <= tol) continue;
    const Real d = segDist(i, i + 1);
    if (d <= upperSeg + tol) continue;

    const std::size_t num =
        static_cast<std::size_t>(std::ceil(d / upperSeg)) - 1;
    if (num == 0) continue;
    const Real dt = (params[i + 1] - params[i]) / (num + 1);
    std::vector<Real> mids;
    mids.reserve(num);
    for (std::size_t k = 1; k <= num; ++k) {
      const Real t = params[i] + dt * static_cast<Real>(k);
      if (t <= params[i] + tol || t >= params[i + 1] - tol) continue;
      mids.push_back(t);
    }
    if (!mids.empty()) pending.push_back({i, std::move(mids)});
  }

  if (pending.empty()) return false;

  std::size_t totalInsert = 0;
  for (const auto &ins : pending) totalInsert += ins.mids.size();

  std::vector<Real> newParams;
  std::vector<int> newFixed;
  newParams.reserve(params.size() + totalInsert);
  newFixed.reserve(fixed.size() + totalInsert);

  std::size_t pendingIdx = 0;
  for (std::size_t i = 0; i < params.size(); ++i) {
    newParams.push_back(params[i]);
    newFixed.push_back(fixed[i]);
    if (pendingIdx < pending.size() && pending[pendingIdx].index == i) {
      const auto &mids = pending[pendingIdx].mids;
      newParams.insert(newParams.end(), mids.begin(), mids.end());
      newFixed.insert(newFixed.end(), mids.size(), false);
      ++pendingIdx;
    }
  }

  params.swap(newParams);
  fixed.swap(newFixed);

  return true;
}

template <int Dim, int Order>
static bool pruneByLower(const Curve<Dim, Order> &curve,
                         const std::vector<Real> &lowerBound,
                         std::vector<Real> &params,
                         std::vector<int> &fixed, Real tol) {
  if (params.size() < 3) return false;

  auto segDist = [&](std::size_t i, std::size_t j) {
    return norm(curve(params[j]) - curve(params[i]), 2);
  };

  std::vector<int> erased;
  for (std::size_t i = 1; i + 1 < params.size(); ++i) {
    if (fixed[i]) {
      continue;
    } 
    const Real prevD = segDist(i - 1, i);
    const Real nextD = segDist(i, i + 1);
    const Real boundPrev =
        std::max<Real>(0, std::min(lowerBound[i - 1], lowerBound[i]));
    const Real boundNext =
        std::max<Real>(0, std::min(lowerBound[i], lowerBound[i + 1]));
    if (prevD < boundPrev - tol || nextD < boundNext - tol) {
      erased.push_back(i);
      ++i;
    }
  }
  if (erased.empty()) return false;
  
  std::vector<Real> newParams;
  std::vector<int> newFixed;
  auto iter = erased.begin();
  for (std::size_t i = 0; i < params.size(); ++i) {
    if (i == *iter) {
      ++iter;
    } else {
      newParams.push_back(params[i]);
      newFixed.push_back(fixed[i]);
    }
  }
  params.swap(newParams);
  fixed.swap(newFixed);

  return true;
}
}  // namespace

template <int Dim, int Order>
OPTNONE_FUNC
auto Curve<Dim, Order>::curvatureBoundParams(
    const std::vector<Real> &params, std::vector<int> &fixed, Real base,
    const std::function<void(const std::vector<Real> &curv,
                             std::vector<Real> &hL, Real &lowerScale,
                             Real &upperScale, Real)> &distBounds) const
    -> std::vector<Real> {
  if constexpr (Order <= 2) {
    throw std::runtime_error(
        "curvatureBoundParams requires curve order greater than 2.");
  }
  if (params.size() < 2) return params;
  std::vector<Real> res = params;
  if (empty() || params.empty() || !distBounds) return res;

  const Real lo = knots.front();
  const Real hi = knots.back();
  const Real tol = distTol();

  // clamp to domain, sort and deduplicate; mark all provided params as fixed
  if (fixed.empty()) {
    std::vector<Real> sorted(params.size());
    std::transform(params.begin(), params.end(), sorted.begin(),
                   [&](Real t) { return std::min(std::max(t, lo), hi); });
    std::sort(sorted.begin(), sorted.end());
    fixed.reserve(sorted.size());
    for (const auto &t : sorted) {
      if (res.empty() || t - res.back() > tol) {
        res.push_back(t);
        fixed.push_back(true);
      }
    }
  }

  auto computeBounds = [&]() {
    std::vector<Real> lower(res.size());
    std::vector<Real> upper(res.size());
    std::vector<Real> hL(res.size());
    Real lowerScale = 0;
    Real upperScale = 0;
    auto curv = curvature(res);
    distBounds(curv, hL, lowerScale, upperScale, base);
    if (hL.size() != res.size()) {
      throw std::runtime_error("distBounds returns size mismatch.");
    }
    lowerScale = 0.2;
    upperScale = 0.6;
    for (std::size_t i = 0; i < res.size(); ++i) {
      lower[i] = std::max<Real>(0, hL[i] * lowerScale);
      upper[i] = hL[i] * upperScale;
    }
    return std::pair<std::vector<Real>, std::vector<Real>>(
        std::move(lower), std::move(upper));
  };

  constexpr std::size_t kMaxIter = 10;

  // enrich by upper distance bound
  std::size_t iter = 0;
  while (iter++ < kMaxIter) {
    auto bounds = computeBounds();
    auto &upperBound = bounds.second;
    if (res.size() < 2) break;

    const bool inserted = insertByUpper(*this, upperBound, res, fixed, tol);
    if (!inserted) break;
  }

  // prune by lower distance bound (never remove initial params)
  iter = 0;
  while (iter++ < kMaxIter) {
    auto bounds = computeBounds();
    auto &lowerBound = bounds.first;
    const bool erased = pruneByLower(*this, lowerBound, res, fixed, tol);
    if (!erased) break;
  }

  return res;
}

template <int Dim, int Order>
Real Curve<Dim, Order>::curvature(const Polynomial<Order, Real> &xPoly,
                                  const Polynomial<Order, Real> &yPoly,
                                  Real t) {
  if (Order <= 2) throw std::runtime_error("Order must be greater than 2.");
  Real xp;
  Real xpp;
  Real yp;
  Real ypp;
  if (t < distTol()) {
    xp = xPoly[1];
    xpp = xPoly[2] * 2;
    yp = yPoly[1];
    ypp = yPoly[2] * 2;
  } else {
    auto xPolyDer = xPoly.der();
    xp = xPolyDer(t);
    xpp = xPolyDer.der()(t);
    auto yPolyDer = yPoly.der();
    yp = yPolyDer(t);
    ypp = yPolyDer.der()(t);
  }
  using std::abs;
  using std::pow;
  using std::sqrt;
  Real numerator = abs(xp * ypp - yp * xpp);
  Real denominator = sqrt(pow((pow(xp, 2) + pow(yp, 2)), 3));
  return numerator / denominator;
}

template <int Order>
OPTNONE_FUNC
void checkFitCurve(const Curve<2, Order> &crv,
                   typename Curve<2, Order>::BCType type,
                   const Vec<Real, 2> &start, const Vec<Real, 2> &end) {
  const auto &knots = crv.getKnots();
  const auto &polys = crv.getPolys();
  auto num = polys.size();
  std::vector<Real> ret(4, 0);
  for (int i = 1; i < num; ++i) {
    ret[0] = std::max(
        ret[0], norm(polys[i - 1](knots[i] - knots[i - 1]) - polys[i][0], 2));
    if constexpr (Order > 2) {
      auto ployDs = der(crv).getPolys();
      auto polyDDs = der(der(crv)).getPolys();
      ret[1] = std::max(
          ret[1],
          norm(ployDs[i - 1](knots[i] - knots[i - 1]) - ployDs[i][0], 2));
      ret[2] = std::max(
          ret[2],
          norm(polyDDs[i - 1](knots[i] - knots[i - 1]) - polyDDs[i][0], 2));
    }
  }
  if constexpr (Order > 2) {
    auto ployDs = der(crv).getPolys();
    auto polyDDs = der(der(crv)).getPolys();
    if (type == Curve<2, 4>::periodic) {
      ret[3] = std::max(
          ret[3],
          norm(polys[num - 1](knots[num] - knots[num - 1]) - polys[0][0], 2));
      ret[3] = std::max(
          ret[3],
          norm(ployDs[num - 1](knots[num] - knots[num - 1]) - ployDs[0][0], 2));
      ret[3] = std::max(
          ret[3],
          norm(polyDDs[num - 1](knots[num] - knots[num - 1]) - polyDDs[0][0],
               2));
    } else if (type == Curve<2, 4>::notAKnot) {
      ret[3] = std::max(ret[3], norm(polys[num - 1][3] - polys[num - 2][3], 2));
      ret[3] = std::max(ret[3], norm(polys[1][3] - polys[0][3], 2));
    }
  }

  std::cout << fmt::v11::format(
      "max error of position: {}, velocity: {}, acceleration: {}, endpoint: "
      "{}\n",
      ret[0], ret[1], ret[2], ret[3]);
}

//=================================================
// explicit instantiation of the followings

template class Curve<2, 1>;
template class Curve<2, 2>;
template class Curve<2, 3>;
template class Curve<2, 4>;
template class Curve<2, 6>;

template Curve<2, 3> der(const Curve<2, 4> &);
template Curve<2, 2> der(const Curve<2, 3> &);
template Curve<2, 1> der(const Curve<2, 2> &);

template Interval<2> boundingBox(const Curve<2, 2> &);
template Interval<2> boundingBox(const Curve<2, 4> &);

template Real area(const Curve<2, 2> &gon, Real tol);
template Real area(const Curve<2, 4> &gon, Real tol);

template Real arclength(const Curve<2, 2> &);
template Real arclength(const Curve<2, 4> &);

template Curve<2, 2> createRect(const Vec<Real, 2> &, const Vec<Real, 2> &);
template Curve<2, 4> createRect(const Vec<Real, 2> &, const Vec<Real, 2> &);

template Curve<2, 2> createLineSegment(const Vec<Real, 2> &p0,
                                       const Vec<Real, 2> &p1);
template Curve<2, 4> createLineSegment(const Vec<Real, 2> &p0,
                                       const Vec<Real, 2> &p1);

template void checkFitCurve<2>(const Curve<2, 2> &crv,
                               typename Curve<2, 2>::BCType type,
                               const Vec<Real, 2> &start,
                               const Vec<Real, 2> &end);
template void checkFitCurve<4>(const Curve<2, 4> &crv,
                               typename Curve<2, 4>::BCType type,
                               const Vec<Real, 2> &start,
                               const Vec<Real, 2> &end);
