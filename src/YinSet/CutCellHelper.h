#include "YinSet/PastingMap.h"
#include "YinSet/PointsLocater.h"
#include "YinSet/YinSet.h"

template <int Order>
struct CutCellHelper {
  using iVec = Vec<int, 2>;
  using rVec = Vec<Real, 2>;
  using Crv = Curve<2, Order>;
  using SRS = SegmentedRealizableSpadjor<Order>::SRS;
  using YinSetPtr = YinSet<2, Order>::YinSetPtr;
  enum CellType { UNKNOWN = -2, OUTER = -1, BOUNDARY = 0, INNER = 1 };
  static auto intersectGridLine(
      const Vec<Real, 2> &lo, const Vec<Real, 2> &hi, const Vec<Real, 2> &h,
      const std::vector<OrientedJordanCurve<2, Order>> &orientedJordanCurves,
      Real tol);

  // partition the curves to vector<vector<<Curve<Dim, Order>>>
  static auto splitCurves(
      rVec lo, rVec h, std::vector<set<Real>> &intersections,
      const std::vector<OrientedJordanCurve<2, Order>> &orientedJordanCurves,
      vector<vector<vector<Curve<2, Order>>>> &gridCurves, Real tol);

  // past curves in every cell.
  static auto pastCells(
      rVec lo, rVec h, Box<DIM> box,
      const vector<vector<vector<Curve<2, Order>>>> &gridCurves,
      vector<vector<YinSetPtr>> &gridSRS, Real tol);

  // fill inner grid.
  static auto fillInner(
      rVec lo, rVec h, Box<DIM> box, const SRS &srs,
      const vector<vector<vector<Curve<2, Order>>>> &gridCurves,
      vector<vector<YinSetPtr>> &gridSRS, Real tol, bool addInner = false);
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
    auto monotonicCrv = jordanCrv.makeMonotonic(tol);
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
        if (std::fabs(p0 - (lo + h * startRow)) < tol) {
          intersectionsForCurve.insert(!swapped ? knotStart : knotsEnd);
          ++startRow;
        }
        if (std::fabs(p1 - (lo + h * endRow)) < tol) {
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
    vector<vector<vector<Curve<2, Order>>>> &gridCurves, Real tol) {
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
      gridCurves[id[0]][id[1]].emplace_back(crv);
    }
  }
}

template <int Order>
// __attribute__((optnone))
auto CutCellHelper<Order>::pastCells(
    rVec lo, rVec h, Box<DIM> box,
    const vector<vector<vector<Curve<2, Order>>>> &gridCurves,
    vector<vector<YinSetPtr>> &gridSRS, Real tol) {
  loop_box_2(box, i0, i1) {
    gridSRS[i0][i1] = nullptr;
    if (gridCurves[i0][i1].empty()) {
      continue;
    }

    auto localLo = lo + h * iVec{i0, i1};
    auto localHi = localLo + h;
    Curve<2, Order> rect = createRect<Order>(localLo, localHi);

    vector<Real> brks;
    Real pre0 = 0;
    Real pre1 = pre0 + h[0];
    Real pre2 = pre1 + h[1];
    Real pre3 = pre2 + h[0];
    PastingMap<Order> pasting(tol * 2);
    for (auto &crv : gridCurves[i0][i1]) {
      auto points = {crv.startpoint(), crv.endpoint()};
      for (auto p : points) {
        if (std::fabs(p[1] - localLo[1]) < tol) {
          p[1] = localLo[1];
          if (std::fabs(p[0] - localLo[0]) < tol) {
            p[0] = localLo[0];
            brks.push_back(pre0);
          } else if (std::fabs(p[0] - localHi[0]) < tol) {
            p[0] = localHi[0];
            brks.push_back(pre1);
          } else {
            brks.push_back(pre0 + p[0] - localLo[0]);
          }
        } else if (std::fabs(p[1] - localHi[1]) < tol) {
          p[1] = localHi[1];
          if (std::fabs(p[0] - localLo[0]) < tol) {
            p[0] = localLo[0];
            brks.push_back(pre3);
          } else if (std::fabs(p[0] - localHi[0]) < tol) {
            p[0] = localHi[0];
            brks.push_back(pre2);
          } else {
            brks.push_back(pre2 + localHi[0] - p[0]);
          }
          brks.push_back(pre2 + localHi[0] - p[0]);
        } else if (std::fabs(p[0] - localLo[0]) < tol) {
          p[0] = localLo[0];
          brks.push_back(pre3 + localHi[1] - p[1]);
        } else if (std::fabs(p[0] - localHi[0]) < tol) {
          p[0] = localHi[0];
          brks.push_back(pre1 + p[1] - localLo[1]);
        }

        pasting.addVertex(p);
      }
    }
    std::sort(brks.begin(), brks.end());
    std::vector<Crv> rectCrvs;
    rect.split(brks, rectCrvs, tol);

    for (auto &crv : gridCurves[i0][i1]) {
      pasting.addEdge(crv, true);
    }
    for (auto &crv : rectCrvs) {
      pasting.addEdge(crv, false);
    }

    vector<OrientedJordanCurve<2, Order>> res;
    pasting.formClosedLoops(res);
    if (res.empty()) continue;
    vector<OrientedJordanCurve<2, Order>> vecCrv;
    for (auto &&crv : res) {
      auto vol = area(crv);
      if (std::abs(vol) > tol) vecCrv.emplace_back(std::move(crv));
    }
    // gridSRS[i0][i1]  = (new YinSet<2, Order>(SRS(vecCrv), tol));
    if (vecCrv.empty()) continue;

    auto yinsetPtr = std::make_shared<YinSet<2, Order>>(SRS(vecCrv), tol);
    if (!yinsetPtr->isBounded()) {
      vecCrv.push_back(rect);
      yinsetPtr =
          std::make_shared<YinSet<2, Order>>(SRS(std::move(vecCrv)), tol);
    }
    gridSRS[i0][i1] = yinsetPtr;
  }
}

template <int Order>
auto CutCellHelper<Order>::fillInner(
    rVec lo, rVec h, Box<DIM> box, const SRS &srs,
    const vector<vector<vector<Curve<2, Order>>>> &gridCurves,
    vector<vector<YinSetPtr>> &gridSRS, Real tol, bool addInner) {
      const int N0 = box.hi()[0] - box.lo()[0] + 1;
  const int N1 = box.hi()[1] - box.lo()[1] + 1;
  vector<vector<int>> tags(N0, vector<int>(N1, CellType::UNKNOWN));
  loop_box_2(box, i0, i1) {
    const auto &boundaryPtr = gridCurves[i0][i1];
    if (!boundaryPtr.empty()) {
      tags[i0][i1] = CellType::BOUNDARY;
    } else {
      tags[i0][i1] = CellType::UNKNOWN;
    }
  }
  auto iLo = box.lo();
  auto iHi = box.hi();
  vector<iVec> queue;
  std::function<void(CellType)> dfsMark = [&queue, &tags, iLo,
                                           iHi](CellType tag) {
    while (!queue.empty()) {
      iVec cur = queue.back();
      queue.pop_back();
      tags[cur[0]][cur[1]] = tag;
      if (cur[0] > iLo[0] && tags[cur[0] - 1][cur[1]] == CellType::UNKNOWN) {
        queue.push_back(iVec{cur[0] - 1, cur[1]});
      }
      if (cur[0] < iHi[0] && tags[cur[0] + 1][cur[1]] == CellType::UNKNOWN) {
        queue.push_back(iVec{cur[0] + 1, cur[1]});
      }
      if (cur[1] > iLo[1] && tags[cur[0]][cur[1] - 1] == CellType::UNKNOWN) {
        queue.push_back(iVec{cur[0], cur[1] - 1});
      }
      if (cur[1] < iHi[1] && tags[cur[0]][cur[1] + 1] == CellType::UNKNOWN) {
        queue.push_back(iVec{cur[0], cur[1] + 1});
      }
    }
  };

  loop_box_2(box, i0, i1) {
    if (tags[i0][i1] == CellType::UNKNOWN) {
      queue.push_back(iVec{i0, i1});
      rVec mid = lo + h * iVec{i0, i1} + h * 0.5;
      int tag = PointsLocater(tol).operator()(srs.getOrientedJordanCurvesRef(),
                                              {mid}, srs.isBounded(tol))[0];
      while (!queue.empty()) dfsMark(CellType(tag));
    }
  }
  if (addInner) {
    loop_box_2(box, i0, i1) {
      if (tags[i0][i1] == CellType::INNER) {
        auto localLo = lo + h * iVec{i0, i1};
        auto localHi = localLo + h;
        Curve<2, Order> rect = createRect<Order>(localLo, localHi);
        vector<OrientedJordanCurve<2, Order>> vecCrv = {rect};
        // gridSRS[i0][i1] = (new YinSet<2, Order>(SRS(vecCrv), tol));
        gridSRS[i0][i1] =
            (std::make_shared<YinSet<2, Order>>(SRS(vecCrv), tol));
      }
    }
  }

  return tags;
}