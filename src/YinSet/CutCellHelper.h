#include "YinSet/PastingMap.h"
#include "YinSet/PointsLocater.h"
#include "YinSet/SegmentedRealizableSpadjor.h"

template <int Order>
struct CutCellHelper {
  using iVec = Vec<int, 2>;
  using rVec = Vec<Real, 2>;
  using Crv = Curve<2, Order>;
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
  static auto fillInner(rVec lo, rVec h, const SRS &srs,
                        Tensor<vector<SRS>, 2> &gridSRS, Real tol);
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
        Real p1 = xPoly(knotsEnd - knotStart);
        bool swapped = false;
        if (p0 > p1) {
          std::swap(p0, p1);
          swapped = true;
        }
        int startRow = std::ceil((p0 - lo - tol) / h);
        int endRow = std::floor((p1 - lo + tol) / h);
        if (std::fabs(p0 - lo + h * startRow) < tol) {
          intersectionsForCurve.insert(!swapped ? knotStart : knotsEnd);
          ++startRow;
        }
        if (std::fabs(p1 - lo + h * endRow) < tol) {
          intersectionsForCurve.insert(!swapped ? knotsEnd : knotStart);
          --endRow;
        }

        for (int r = startRow; r <= endRow; ++r) {
          Real row = lo + r * h;
          auto intersection =
              root(xPoly - row, (knotsEnd - knotStart) / 2, tol, newtonMaxIter);
          if (intersection > 0 && intersection < knotsEnd - knotStart)
            intersectionsForCurve.insert(intersection + knotStart);
        }
      };

      for (int k = 0; k < DIM; ++k) {
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
        if (vcmp(prevPoint, polys[j](*it - knots[j])) == 0) {
          it = intersectionsForCurve.erase(it);
        } else {
          prevPoint = polys[j](*it - knots[j]);
          ++it;
        }
      }
      --it;
      prevPoint = polys[(j + 1) % numPolys][0];
      if (vcmp(prevPoint, polys[j](*it - knots[j])) == 0) {
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
    // auto &polys = crv.getPolys();
    // auto &knots = crv.getKnots();

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
  loop_box_2(gridCurves.box(), i0, i1) {
    auto localLo = lo + h * iVec{i0, i1};
    auto localHi = localLo + h;
    Curve<2, Order> rect = createRect<Order>(localLo, localHi);

    vector<Real> brks;
    Real pre0 = 0;
    Real pre1 = pre0 + h[0];
    Real pre2 = pre1 + h[1];
    Real pre3 = pre2 + h[0];
    for (auto &crv : gridCurves(i0, i1)) {
      auto points = {crv.startpoint(), crv.endpoint()};
      for (auto &p : points) {
        if (std::fabs(p[1] - localLo[1]) < tol)
          brks.push_back(pre0 + p[0] - localLo[0]);
        else if (std::fabs(p[1] - localHi[1]) < tol)
          brks.push_back(pre2 + localHi[0] - p[0]);
        else if (std::fabs(p[0] - localLo[0]) < tol)
          brks.push_back(pre3 + localHi[1] - p[1]);
        else if (std::fabs(p[0] - localHi[0]) < tol)
          brks.push_back(pre1 + p[1] - localLo[1]);
      }
    }
    std::sort(brks.begin(), brks.end());
    std::vector<Crv> rectCrvs;
    rect.split(brks, rectCrvs, tol);

    PastingMap<Order> pasting(tol);
    for (auto &crv : gridCurves(i0, i1)) {
      pasting.addEdge(crv, true);
    }
    for (auto &crv : rectCrvs) {
      pasting.addEdge(crv, false);
    }

    vector<Curve<2, Order>> res;
    pasting.formClosedLoops(res);
    for (auto& crv : res) {
      gridSRS(i0, i1).emplace_back(vector<OrientedJordanCurve<2, Order>>{crv});
    }
  }
}

template <int Order>
auto CutCellHelper<Order>::fillInner(rVec lo, rVec h, const SRS &srs,
                                     Tensor<vector<SRS>, 2> &gridSRS,
                                     Real tol) {
  Tensor<int, 2> tags(gridSRS.box());
  enum { UNKNOWN = -2, OUTER = -1, BOUNDARY = 0, INNER = 1 };
  loop_box_2(gridSRS.box(), i0, i1) {
    auto &srs = gridSRS(i0, i1);
    if (srs.empty()) {
      tags(i0, i1) = UNKNOWN;
    } else {
      tags(i0, i1) = BOUNDARY;
    }
  }
  auto iLo = gridSRS.box().lo();
  auto iHi = gridSRS.box().hi();
  vector<iVec> queue;
  std::function<void(int)> dfsMark = [&queue, &tags, iLo, iHi](int tag) {
    while (!queue.empty()) {
      iVec cur = queue.back();
      queue.pop_back();
      tags(cur) = tag;
      if (cur[0] > iLo[0] && tags(cur[0] - 1, cur[1]) == UNKNOWN) {
        queue.push_back(iVec{cur[0] - 1, cur[1]});
      }
      if (cur[0] < iHi[0] && tags(cur[0] + 1, cur[1]) == UNKNOWN) {
        queue.push_back(iVec{cur[0] + 1, cur[1]});
      }
      if (cur[1] > iLo[1] && tags(cur[0], cur[1] - 1) == UNKNOWN) {
        queue.push_back(iVec{cur[0], cur[1] - 1});
      }
      if (cur[1] < iHi[1] && tags(cur[0], cur[1] + 1) == UNKNOWN) {
        queue.push_back(iVec{cur[0], cur[1] + 1});
      }
    }
  };

  loop_box_2(tags.box(), i0, i1) {
    if (tags(i0, i1) == UNKNOWN) {
      queue.push_back(iVec{i0, i1});
      rVec mid = lo + h * iVec{i0, i1} + h * 0.5;
      int tag = PointsLocater(tol).operator()(srs.getOrientedJordanCurvesRef(),
                                              {mid})[0];
      while (!queue.empty()) dfsMark(tag);
    }
  }

  loop_box_2(tags.box(), i0, i1) {
    if (tags(i0, i1) == INNER) {
      auto localLo = lo + h * iVec{i0, i1};
      auto localHi = localLo + h;
      Curve<2, Order> rect = createRect<Order>(localLo, localHi);
      gridSRS(i0, i1).emplace_back(
          std::vector<OrientedJordanCurve<2, Order>>{rect});
    }
  }
}