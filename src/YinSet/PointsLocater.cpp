#include "PointsLocater.h"

#include <algorithm>
#include <limits>

auto PointsLocater::compute(const vector<Segment<2>> &segs,
                            const vector<rVec> &pts, int defVal)
    -> vector<int> {
  vector<int> output(pts.size());

  // reserve space for the back references
  backref.resize(segs.size());

  // initialize the event queue
  // for those endpoints
  for (auto sit = segs.begin(); sit != segs.end(); ++sit) {
    int type = pt_cmp.compare(sit->p[0], sit->p[1]);
    assert(type != 0);
    event.emplace_back(sit->p[0], type, -1, sit);
    event.emplace_back(sit->p[1], -type, -1, sit);
  }
  // for the query points
  for (auto pit = pts.begin(); pit != pts.end(); ++pit)
    event.emplace_back(*pit, -3, pit - pts.begin(), segs.end());

  // minor keyword : lower < upper < query
  auto sorter = [this](const Event &e1, const Event &e2) -> bool {
    int r = pt_cmp.compare(e1.p, e2.p);
    return (r == 0) ? (e1.type > e2.type) : (r == -1);
  };
  std::sort(event.begin(), event.end(), sorter);

  int rem = pts.size();
  for (int n = 0; rem > 0; n++) {
    const Event &evt = event[n];
    // inform the compare object
    ep = evt.p;

    if (evt.type == -1) {  // upper
      auto r = status.insert(evt.cit);
      assert(r.second);
      backref[evt.cit - segs.begin()] = r.first;
    } else if (evt.type == 1) {  // lower
      status.erase(backref[evt.cit - segs.begin()]);
    } else {  // query
      // deal with the case that query point coincides with endpoint
      if (n >= 1 && norm(event[n].p - event[n - 1].p) < tol) {
        if (event[n - 1].type == -3)  // also a query point ?
          output[event[n].outIdx] = output[event[n - 1].outIdx];
        else
          output[event[n].outIdx] = 0;
        --rem;
        continue;
      }

      // construct a dummy segment vertical segment
      rVec q = evt.p;
      vector<Segment<2>> dummy;
      dummy.emplace_back(q, rVec{q[0], q[1] - 1e3 * tol});
      auto r = status.insert(dummy.begin());

      if (!r.second) {
        output[evt.outIdx] = 0;
      } else {
        // grab a neighbouring segment
        StatusType::const_iterator i = r.first;
        StatusType::const_iterator left = status.end(), right = status.end();

        do {
          // check the left neighbour
          if (i != status.begin()) {
            left = i;
            --left;
            if ((*left)->contain(q, tol)) {
              output[evt.outIdx] = 0;
              break;
            }
          }
          // check the right neighbour
          right = i;
          ++right;
          if (right != status.end()) {
            if ((*right)->contain(q, tol)) {
              output[evt.outIdx] = 0;
              break;
            }
          }

          // no neighbouring segment
          if (left == status.end() && right == status.end()) {
            output[evt.outIdx] = defVal;  // use default value
            break;
          }
          // a neighbouring segment is found and not ON
          if (right == status.end()) {
            output[evt.outIdx] = ((*left)->area(q) > 0) ? (1) : (-1);
          } else {
            output[evt.outIdx] = ((*right)->area(q) > 0) ? (1) : (-1);
          }
        } while (false);

        status.erase(i);
      }  // end insertion success
      --rem;
    }  // end query
  }    // end for

  // early exit is highly possible, so free the resource
  backref.clear();
  event.clear();
  status.clear();
  return output;
}

template <int Order>
vector<int> PointsLocater::operator()(
    const vector<OrientedJordanCurve<DIM, Order>> &ys,
    const vector<rVec> &queries, bool bounded) {
  std::vector<int> ret;
  if (ys.empty() || queries.empty()) {
    return {};
  }
  using Crv = Curve<DIM, Order>;

  // <distance, point, normal, direction>
  std::vector<std::tuple<Real, rVec, rVec, Crv>> closest_points;
  closest_points.resize(queries.size(),
                        std::make_tuple(std::numeric_limits<Real>::max(),
                                        rVec{0, 0}, rVec{0, 0}, Crv()));
  const int ix = 0;
  const int iy = 1;
  VecCompare<Real, 2> vcmp(tol);

  for (const auto &crv : ys) {
    auto monotonicCurve = crv.makeMonotonic(tol);
    auto &polys = monotonicCurve.getPolys();
    auto &knots = monotonicCurve.getKnots();
    for (int i = 0; i < polys.size(); ++i) {
      auto &poly = polys[i];
      auto xPoly = getComp(poly, ix);
      auto p0 = poly[0];
      auto p1 = poly(knots[i + 1] - knots[i]);
      Crv localCrv(poly, knots[i + 1] - knots[i]);
      if (p0[ix] > p1[ix]) std::swap(p0, p1);

      for (int j = 0; j < queries.size(); ++j) {
        const auto &q = queries[j];
        Real x = q[ix];
        if (x <= p0[ix] - tol || x >= p1[ix] + tol) continue;
        if (norm(p0[iy] - q[iy]) > std::get<0>(closest_points[j]) &&
            norm(p1[iy] - q[iy]) > std::get<0>(closest_points[j]))
          continue;

        auto brk =
            root(xPoly - x, (knots[i + 1] - knots[i]) / 2, tol, newtonMaxIter);
        if (brk < tol || brk > knots[i + 1] - knots[i] - tol) {
          if (brk < tol && brk > -tol)
            brk = 0;
          else if (brk > knots[i + 1] - knots[i] - tol &&
                   brk < knots[i + 1] - knots[i] + tol)
            brk = knots[i + 1] - knots[i];
          else
            throw std::runtime_error("root not found");
        }
        auto y = poly(brk)[iy];
        auto d = norm(y - q[iy]);
        rVec aim = {x, y};
        if (norm(d - std::get<0>(closest_points[j])) < tol) {
          auto preCrv = std::get<3>(closest_points[j]);
          auto compareRes = localCrv.compare(preCrv, q, ix, iy, tol);
          if (compareRes * (y - q[iy]) > 0) {
            closest_points[j] = std::make_tuple(
                d, aim, localCrv.normalVector(brk, tol)[0], localCrv);
          }
        } else if (d < std::get<0>(closest_points[j])) {
          closest_points[j] = std::make_tuple(
              d, aim, localCrv.normalVector(brk, tol)[0], localCrv);
        }
      }
    }
  }

  ret.resize(queries.size());
  for (int i = 0; i < queries.size(); ++i) {
    auto &cp = closest_points[i];
    if (std::get<0>(cp) == std::numeric_limits<Real>::max()) {
      ret[i] = bounded ? -1 : 1;
    } else if (std::get<0>(cp) > tol) {
      ret[i] = dot(std::get<2>(cp), std::get<1>(cp) - queries[i]) > 0 ? 1 : -1;
    } else {
      ret[i] = 0;
    }
  }
  return ret;
}

template vector<int> PointsLocater::operator()<4>(
    const vector<OrientedJordanCurve<2, 4>> &ys, const vector<rVec> &queries,
    bool);
template vector<int> PointsLocater::operator()<2>(
    const vector<OrientedJordanCurve<2, 2>> &ys, const vector<rVec> &queries,
    bool);