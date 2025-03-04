#include "YinSet/SegmentedRealizableSpadjor.h"

template <int Order>
struct CutCellHelper {
  using rVec = Vec<Real, 2>;
  using SRS = SegmentedRealizableSpadjor<Order>::SRS;
  static auto intersectGridLine(
      const Vec<Real, 2> &lo, const Vec<Real, 2> &hi, const Vec<Real, 2> &h,
      const std::vector<OrientedJordanCurve<2, Order>> &orientedJordanCurves,
      Real tol);

  // partition the curves to Tensor<Curve<Dim, Order>, 2>
  static auto splitCurves(
      rVec lo, rVec h, std::vector<set<Real>> &intersections,
      const std::vector<OrientedJordanCurve<2, Order>> &orientedJordanCurves,
      Tensor<vector<Curve<2, Order>>, 2> &gridCurves, Real tol);

  // past curves in every cell.
  static auto pastCells(rVec lo, rVec h,
                        const Tensor<vector<Curve<2, Order>>, 2> &gridCurves,
                        Tensor<vector<SRS>, 2> &gridSRS, Real tol);

  // fill inner grid.
  static auto fillInner(const SRS& srs, Tensor<vector<SRS>, 2> &gridSRS, Real tol);
};

template <int Order>
auto CutCellHelper<Order>::intersectGridLine(
    const Vec<Real, 2> &lo, const Vec<Real, 2> &hi, const Vec<Real, 2> &h,
    const std::vector<OrientedJordanCurve<2, Order>> &orientedJordanCurves,
    Real tol) {
  using std::vector;
  vector<set<Real>> intersections;

  for (const auto &jordanCrv : orientedJordanCurves) {
    intersections.emplace_back();
    auto &intersectionsForCurve = intersections.back();
    auto monotonicCrv = jordanCrv.makeMonotonic(distTol());
    auto knots = monotonicCrv.getKnots();
    auto polys = monotonicCrv.getPolys();
    auto numPolys = polys.size();
    for (int i = 0; i < numPolys; ++i) {
      // lambda for xRow and yCol
      auto intersect = [tol](int k, Real lo, const auto &xPoly, Real h,
                             Real knotStart, Real knotsEnd,
                             auto &intersectionsForCurve, int newtonMaxIter) {
        Real p0 = xPoly[0];
        Real p1 = xPoly(knotsEnd);
        bool swapped = false;
        if (p0 > p1) {
          std::swap(p0, p1);
          swapped = true;
        }
        int startRow = std::ceil((p0 - lo - tol) / h);
        int endRow = std::floor((p1 - lo + tol) / h);
        if (std::abs(p0 - lo + h * startRow) < tol) {
          intersectionsForCurve.insert(!swapped ? knotStart : knotsEnd);
          ++startRow;
        }
        if (std::abs(p1 - lo + h * endRow) < tol) {
          intersectionsForCurve.insert(!swapped ? knotsEnd : knotStart);
          --endRow;
        }

        for (int r = startRow; r <= endRow; ++r) {
          Real row = lo + r * h;
          auto intersection = root(xPoly - row, 0, tol, newtonMaxIter);
          intersectionsForCurve.insert(intersection);
        }
      };

      for (int k = 0; k < 2; ++k) {
        RealPolynomial<Order> xPoly = getComp(polys[i], k);
        intersect(k, lo[k], xPoly, h[k], knots[i], knots[i + 1],
                  intersectionsForCurve, newtonMaxIter);
      }
    }
  }

  return intersections;
}

template <int Order>
auto CutCellHelper<Order>::splitCurves(
    rVec lo, rVec h, std::vector<set<Real>> &intersections,
    const std::vector<OrientedJordanCurve<2, Order>> &orientedJordanCurves,
    Tensor<vector<Curve<2, Order>>, 2> &gridCurves, Real tol) {
  using Point = Vec<Real, 2>;
  VecCompare<Real, 2> vcmp(tol);
  // remove too close intersections
  for (int i = 0; i < intersections.size(); ++i) {
    auto &intersectionsForCurve = intersections[i];
    auto &crv = orientedJordanCurves[i];
    auto &polys = crv.getPolys();
    auto &knots = crv.getKnots();
    size_t numPolys = polys.size();
    Point prevPoint = polys[0][0];
    auto it = intersectionsForCurve.begin();
    for (int j = 0; j < numPolys; ++j) {
      while (it != intersectionsForCurve.end() && *it < knots[j + 1]) {
        if (vcmp(prevPoint, polys[j](*it)) == 0) {
          it = intersectionsForCurve.erase(it);
        } else {
          prevPoint = polys[j](*it);
          ++it;
        }
      }
      --it;
      prevPoint = polys[(j + 1) % numPolys][0];
      if (vcmp(prevPoint, polys[j](*it)) == 0) {
        it = intersectionsForCurve.erase(it);
      } else {
        ++it;
      }
    }
  }

  // construct split curve on intersections
  for (int i = 0; i < intersections.size(); ++i) {
    auto &intersectionsForCurve = intersections[i];
    auto &crv = orientedJordanCurves[i];
    auto &polys = crv.getPolys();
    auto &knots = crv.getKnots();
    size_t numPolys = polys.size();

    vector<Curve<2, Order>> splitRes;
    vector<Real> splitKnots(intersectionsForCurve.begin(),
                            intersectionsForCurve.end());
    crv.split(splitKnots, splitRes, tol);

    for (auto &&crv : splitRes) {
      rVec mid = crv.midpoint();
      auto id = Vec<int, 2>((mid - lo) / h);
      gridCurves(id).emplace_back(std::move(crv));
    }
  }
}

template <int Order>
auto CutCellHelper<Order>::pastCells(
    rVec lo, rVec h, const Tensor<vector<Curve<2, Order>>, 2> &gridCurves,
    Tensor<vector<SRS>, 2> &gridSRS, Real tol) {

    }

template <int Order>
auto CutCellHelper<Order>::fillInner(const SRS& srs, Tensor<vector<SRS>, 2> &gridSRS, Real tol) {
  Tensor<int, 2> tags(gridSRS.box());
  enum {UNKNOWN = -1, INNER = 0, OUTER = 1, BOUNDARY = 2};
  loop_box_2(gridSRS.box(), i0, i1) {
    auto &srs = gridSRS(i0, i1);
    if (srs.empty()) {
      tags(i0, i1) = UNKNOWN;
    } else {
      tags(i0, i1) = BOUNDARY;
    }
  }

}