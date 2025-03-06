#include "PastingMap.h"

#include <algorithm>

template <int Order, class Selector>
void PastingMap<Order, Selector>::removeEdge(const rVec &oldtail,
                                             std::set<size_t> &repo, Iter eit) {
  necessaryEdge.erase(*eit);
  repo.erase(eit);
  if (repo.empty()) graph.erase(oldtail);
}

template <int Order, class Selector>
void PastingMap<Order, Selector>::formClosedLoops(vector<Crv> &outCont) {
  Crv jordan;
  rVec oldtail;
  rVec newtail;
  vector<std::pair<rVec, Real>> footprint;
  VecCompare<Real, SpaceDim> vcmp(tol);
  // for (auto& crv : allCrvs) {
  //   std::cout << crv.startpoint() << ", " << crv.endpoint() << std::endl;
  // }

  while (!necessaryEdge.empty()) {
    Iter outEdge;
    std::set<size_t> *cands;
    if (jordan.empty()) {  // start a new loop
      auto id = *necessaryEdge.begin();
      auto startpoint = allCrvs[id].startpoint();
      auto iter = graph.find(startpoint);
      cands = &(iter->second);
      outEdge = cands->find(id);
      // take care of the starting point
      oldtail = iter->first;
      footprint.resize(1);
      footprint.front() = std::make_pair(oldtail, jordan.domain().lo()[0]);
    } else {  // grow
      oldtail = newtail;
      assert(graph.find(oldtail) != graph.end());
      cands = &(graph[oldtail]);
      outEdge = std::min_element(cands->cbegin(), cands->cend(),
                                 Selector(tol, oldtail, jordan, allCrvs));
    }
    jordan.concat(allCrvs[*outEdge]);
    newtail = jordan.endpoint();
    // delete the already-used out-edge
    removeEdge(oldtail, *cands, outEdge);

    auto h = std::find_if(
        footprint.cbegin(), footprint.cend(),
        [&](const auto &obj) { return vcmp.compare(newtail, obj.first) == 0; });
    if (h != footprint.cend()) {  // a cycle is formed
      outCont.push_back(
          jordan.extract(h->second, jordan.domain().hi()[0], tol));
      jordan = jordan.extract(jordan.domain().lo()[0], h->second, tol);
      footprint.erase(++h, footprint.cend());
    } else {
      footprint.emplace_back(newtail, jordan.domain().hi()[0]);
    }
  }
  assert(jordan.empty());
}

template <int Order>
bool OutEdgeSelectorByKnots<Order>::operator()(const size_t &lhsId,
                                               const size_t &rhsId) {
  // TODO(ytan): there is a bug in selecting multi curves.
  const auto &lhs = allCrvs[lhsId];
  const auto &rhs = allCrvs[rhsId];
  auto getNextPoint = [](const Crv &p) {
    const auto &knots = p.getKnots();
    return p.getPolys().front()(knots[1] - knots[0]);
  };
  auto lhsPoint = getNextPoint(lhs);
  auto rhsPoint = getNextPoint(rhs);
  // auto [lhsPoint, rhsPoint] = lhs.getComparablePoint(rhs, tol, 0);
  rVec d1 = normalize(lhsPoint - standpoint);
  rVec d2 = normalize(rhsPoint - standpoint);
  Real s1 = cross(indir, d1);
  Real s2 = cross(indir, d2);
  if (s1 * s2 <= 0) return (s1 >= 0);
  Real c1 = dot(indir, d1);
  Real c2 = dot(indir, d2);
  if (s1 > 0) return c1 < c2;
  return c1 > c2;
}

//==========================================================
template class PastingMap<2>;
template class PastingMap<4>;
