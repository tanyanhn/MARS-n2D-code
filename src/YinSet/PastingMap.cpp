#include "PastingMap.h"

#include <algorithm>



template <int Order, class Selector>
void PastingMap<Order, Selector>::removeEdge(const rVec &oldtail,
                                             std::set<size_t> &repo, Iter eit) {
  assert(repo.contains(*eit) && "eit is not in repo");
  necessaryEdge.erase(*eit);
  repo.erase(eit);
  if (repo.empty()) graph.erase(oldtail);
}

template <int Order, class Selector> 
// __attribute__((optnone))
void PastingMap<Order, Selector>::formClosedLoops(vector<OrientedJordanCurve<DIM, Order>> &outCont) {
  Crv jordan;
  rVec oldtail;
  rVec newtail;
  vector<std::pair<rVec, Real>> footprint;
  VecCompare<Real, SpaceDim> vcmp(tol);
  // for (auto& crv : allCrvs) {
  //   std::cout << crv.startpoint() << ", " << crv.endpoint() << std::endl;
  // }

  while (!necessaryEdge.empty() || !jordan.empty()) {
    Iter outEdge;
    size_t tmp;
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
      auto iter = graph.find(oldtail);
      if(iter == graph.end()) breakpoint();
      iter = graph.find(oldtail);
      cands = &(iter->second);
      outEdge = std::min_element(cands->cbegin(), cands->cend(),
                                 Selector(tol, jordan, allCrvs));
      tmp = *outEdge;
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
                                               const size_t &rhsId) const {
  return compare(lhsId, rhsId) == -1;
}
template <int Order>
// __attribute__((optnone))
int OutEdgeSelectorByKnots<Order>::compare(const size_t &lhsId,
                                           const size_t &rhsId) const {
  const auto &lhs = allCrvs[lhsId];
  const auto &rhs = allCrvs[rhsId];
  auto d1 = lhs.getComparableDirection(tol, 0);
  auto d2 = rhs.getComparableDirection(tol, 0);
  if (norm(d1 - d2) < tol) return 0;
  Real s1 = cross(d1, indir);
  Real s2 = cross(d2, indir);
  Real c1 = dot(d1, indir);
  Real c2 = dot(d2, indir);
  Real angle1 = std::atan2(s1, c1);
  if (angle1 < 0) angle1 += 2 * M_PI;
  Real angle2 = std::atan2(s2, c2);
  if (angle2 < 0) angle2 += 2 * M_PI;
  return angle1 <= angle2 ? -1 : 1;
}

//==========================================================
template class PastingMap<2>;
template class PastingMap<4>;
