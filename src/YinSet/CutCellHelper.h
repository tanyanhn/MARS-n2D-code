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
  enum CellType {
    UNKNOWN = -3,
    OUTER = -1,
    BOUNDARY = 0,
    INNER = 1,
    BOUNDOUTER = -2,  // for classify include boundary bug total outer cell.
    BOUNDINNER = 2    // same as above.
  };
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
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
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
    using localReal = long double;
    for (int i = 0; i < numPolys; ++i) {
      // lambda for xRow and yCol
      auto intersect = [tol](int k, localReal lo, const auto &xPoly,
                             localReal h, localReal knotStart,
                             localReal knotsEnd, auto &intersectionsForCurve,
                             int newtonMaxIter)
#ifdef OPTNONE
          __attribute__((optnone))
#endif  // OPTNONE
      {
        localReal p0 = xPoly[0];
        localReal p1 = xPoly(knotsEnd - knotStart);
        bool swapped = false;
        if (p0 > p1) {
          std::swap(p0, p1);
          swapped = true;
        }
        int startRow = std::ceil((p0 - lo - tol) / h);
        int endRow = std::floor((p1 - lo + tol) / h);
        if (std::fabs(p0 - (lo + h * startRow)) <= tol) {
          auto value = !swapped ? knotStart : knotsEnd;
          intersectionsForCurve.insert(value);
          // auto p = xPoly(value);
          ++startRow;
        }
        if (std::fabs(p1 - (lo + h * endRow)) <= tol) {
          auto value = !swapped ? knotsEnd : knotStart;
          intersectionsForCurve.insert(value);
          // auto p = xPoly(value);
          --endRow;
        }

        for (int r = startRow; r <= endRow; ++r) {
          localReal row = lo + r * h;
          auto intersection =
              root(xPoly - row, (knotsEnd - knotStart) / 2, tol, newtonMaxIter);
          if (intersection > 0 && intersection < knotsEnd - knotStart) {
            intersectionsForCurve.insert(intersection + knotStart);
            // auto p = xPoly(intersection);
            // auto k = (p - row);
          }
        }
      };

      for (int k = 0; k < DIM; ++k) {
        Polynomial<Order, localReal> xPoly(getComp(polys[i], k));
        intersect(k, lo[k], xPoly, h[k], knots[i], knots[i + 1],
                  intersectionsForCurve, newtonMaxIter());
      }
    }
  }

  return intersections;
}

template <int Order>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
auto CutCellHelper<Order>::splitCurves(
    rVec lo, rVec h, std::vector<set<Real>> &intersections,
    const std::vector<OrientedJordanCurve<2, Order>> &orientedJordanCurves,
    vector<vector<vector<Curve<2, Order>>>> &gridCurves, Real tol) {
  // construct split curve on intersections
  for (int i = 0; i < intersections.size(); ++i) {
    auto &intersectionsForCurve = intersections[i];
    auto &crv = orientedJordanCurves[i];
    if (!crv.isClosed(distTol()))
      throw std::runtime_error("Curve must be closed.");
    // auto &polys = crv.getPolys();
    // auto &knots = crv.getKnots();

    vector<Curve<2, Order>> splitRes;
    vector<Real> splitKnots(intersectionsForCurve.begin(),
                            intersectionsForCurve.end());
    crv.split(splitKnots, splitRes, tol, true);

    for (auto &&crv : splitRes) {
      // if (crv.empty()) continue;
      // auto start = crv.startpoint();
      // auto end = crv.endpoint();
      rVec mid = crv.midpoint();
      auto id = Vec<int, 2>((mid - lo) / h);
      // auto i0 = id[0];
      // auto i1 = id[1];
      gridCurves[id[0]][id[1]].emplace_back(crv);
      // for (auto v : splitKnots) {
      //   auto p = crv(v);
      //   auto id = Vec<int, 2>((p - lo) / h);
      // }
    }
  }
}

template <int Order>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
auto CutCellHelper<Order>::pastCells(
    rVec lo, rVec h, Box<DIM> box,
    const vector<vector<vector<Curve<2, Order>>>> &gridCurves,
    vector<vector<YinSetPtr>> &gridSRS, Real tol) {
  loop_box_2(box, i0, i1) {
    gridSRS[i0][i1] = nullptr;
    if (gridCurves[i0][i1].empty()) {
      continue;
    }

    using localReal = Real;

    Vec<localReal, 2> localLo(lo + h * iVec{i0, i1});
    Vec<localReal, 2> localHi = localLo + h;
    Curve<2, Order> rect = createRect<Order>(localLo, localHi);

    vector<localReal> brks;
    localReal pre0 = 0;
    localReal pre1 = pre0 + h[0];
    localReal pre2 = pre1 + h[1];
    localReal pre3 = pre2 + h[0];
    // tol * 2 since we move point while attach cell,
    PastingMap<Order, localReal> pasting(tol * 2);
    pasting.setPeriod(pre3 + h[1]);
    vector<OrientedJordanCurve<2, Order>> res;
    for (auto &crv : gridCurves[i0][i1]) {
      if (crv.isClosed(tol)) {
        res.emplace_back(crv);
        continue;
      }
      Polynomial<Order, Vec<long double, DIM>> tmppoly(crv.getPolys()[22]);
      auto tmpT = crv.getKnots()[23] - crv.getKnots()[22];
      auto tmpP = tmppoly(tmpT);
      auto points = {crv.startpoint(), crv.endpoint()};
      for (auto p : points) {
        if (std::fabs(p[1] - localLo[1]) < tol) {
          if (std::fabs(p[0] - localLo[0]) < tol) {
            brks.push_back(pre0);
          } else if (std::fabs(p[0] - localHi[0]) < tol) {
            brks.push_back(pre1);
          } else {
            brks.push_back(pre0 + p[0] - localLo[0]);
          }
        } else if (std::fabs(p[1] - localHi[1]) < tol) {
          if (std::fabs(p[0] - localLo[0]) < tol) {
            brks.push_back(pre3);
          } else if (std::fabs(p[0] - localHi[0]) < tol) {
            brks.push_back(pre2);
          } else {
            brks.push_back(pre2 + localHi[0] - p[0]);
          }
        } else if (std::fabs(p[0] - localLo[0]) < tol) {
          brks.push_back(pre3 + localHi[1] - p[1]);
        } else if (std::fabs(p[0] - localHi[0]) < tol) {
          brks.push_back(pre1 + p[1] - localLo[1]);
        } else {
          throw std::runtime_error("endpoint can't attach to cell.");
        }
      }
      pasting.addCellEdge(std::move(crv), *(++brks.rbegin()), brks.back(), true);
    }
    std::sort(brks.begin(), brks.end());
    std::vector<Crv> rectCrvs;
    rect.split(brks, rectCrvs, tol, true);
    for (auto &&crv : rectCrvs) {
      if (!crv.empty())
        pasting.addCellEdge(std::move(crv), crv.getKnots()[0], crv.getKnots().back(),
                            false);
    }

    pasting.formClosedLoops(res);
    if (res.empty()) continue;
    vector<OrientedJordanCurve<2, Order>> vecCrv;
    for (auto &&crv : res) {
      auto vol = area(crv);
      // TODO(ytan) newtonTol too small add extra miner holes, and add cell
      // inner. distTol too big to ignore small area.
      if (vol > newtonTol() || vol < -newtonTol())
        vecCrv.emplace_back(std::move(crv));
    }
    if (vecCrv.empty()) continue;
    // tol * 2 same as pasting(tol * 2).
    auto yinsetPtr = std::make_shared<YinSet<2, Order>>(SRS(vecCrv), distTol());
    // if (!yinsetPtr->isBounded()) {
    //   vecCrv.push_back(rect);
    //   yinsetPtr =
    //       std::make_shared<YinSet<2, Order>>(SRS(std::move(vecCrv)),
    //       distTol());
    // }
    gridSRS[i0][i1] = yinsetPtr;
  }
}

template <int Order>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
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

  bool bounded = srs.isBounded(0);
  loop_box_2(box, i0, i1) {
    if (tags[i0][i1] == CellType::UNKNOWN) {
      queue.push_back(iVec{i0, i1});
      rVec mid = lo + h * iVec{i0, i1} + h * 0.5;
      int tag = PointsLocater(tol).operator()(srs.getOrientedJordanCurvesRef(),
                                              {mid}, bounded)[0];
      while (!queue.empty()) dfsMark(CellType(tag));
    }
  }
  loop_box_2(box, i0, i1) {
    if (tags[i0][i1] == CellType::BOUNDARY && gridSRS[i0][i1] == nullptr) {
      rVec mid = lo + h * iVec{i0, i1} + h * 0.5;
      int tag = PointsLocater(tol).operator()(srs.getOrientedJordanCurvesRef(),
                                              {mid}, bounded)[0];
      if (tag == CellType::INNER) {
        tags[i0][i1] = CellType::BOUNDINNER;
      } else if (tag == CellType::OUTER) {
        tags[i0][i1] = CellType::BOUNDOUTER;
      } else {
        // may be result of makeMonotonic(to)::translate(), have to enlarge
        // distTol.
        throw std::runtime_error("Boundary cell must be inner or outer.");
      }
    }
    auto localLo = lo + h * iVec{i0, i1};
    auto localHi = localLo + h;
    Curve<2, Order> rect = createRect<Order>(localLo, localHi);
    if (gridSRS[i0][i1] != nullptr) {
      auto vol = gridSRS[i0][i1]->area();
      bool addRect = vol < -1e-8;
      if (std::abs(vol) < distTol()) {
        rVec mid = lo + h * iVec{i0, i1} + h * 0.5;
        int tag = PointsLocater(tol).operator()(
            srs.getOrientedJordanCurvesRef(), {mid}, bounded)[0];
        addRect = (tag == CellType::INNER);
      }
      if (addRect) {
        auto vecCrv = gridSRS[i0][i1]->getBoundaryCycles();
        vecCrv.push_back(rect);
        gridSRS[i0][i1] = std::make_shared<YinSet<2, Order>>(
            SRS(std::move(vecCrv)), distTol());
      }
    }
  }
  if (addInner) {
    loop_box_2(box, i0, i1) {
      if (tags[i0][i1] == CellType::INNER ||
          tags[i0][i1] == CellType::BOUNDINNER) {
        auto localLo = lo + h * iVec{i0, i1};
        auto localHi = localLo + h;
        Curve<2, Order> rect = createRect<Order>(localLo, localHi);
        vector<OrientedJordanCurve<2, Order>> vecCrv = {rect};
        gridSRS[i0][i1] =
            (std::make_shared<YinSet<2, Order>>(SRS(vecCrv), distTol()));
      }
    }
  }

  return tags;
}