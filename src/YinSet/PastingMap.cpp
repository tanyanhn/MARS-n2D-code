#include "PastingMap.h"

#include <algorithm>

template <int Order, typename Num, class Selector>
void PastingMap<Order, Num, Selector>::removeEdge(const rVec &oldtail, auto &repo,
                                             Iter eit) {
  necessaryEdge.erase(*eit);
  repo->second.erase(eit);
  if (repo->second.empty()) graph.erase(oldtail);
}

template <int Order, typename Num, class Selector>
OPTNONE_FUNC
void PastingMap<Order, Num, Selector>::formClosedLoops(vector<OrientedJordanCurve<DIM, Order>> &outCont) {
  if (useCellGraph) return formCellClosedLoops(outCont);
  Crv jordan;
  rVec oldtail;
  rVec newtail;
  vector<std::pair<rVec, Num>> footprint;
  VecCompare<Num, SpaceDim> vcmp(tol);

  while (!necessaryEdge.empty() || !jordan.empty()) {
    Iter outEdge;
    // size_t tmp;
    auto cands = graph.begin();
    if (jordan.empty()) {  // start a new loop
      auto id = *necessaryEdge.begin();
      auto startpoint = allCrvs[id].startpoint();
      auto iter = graph.find(startpoint);
      cands = iter;
      outEdge = cands->second.find(id);
      // take care of the starting point
      oldtail = iter->first;
      footprint.resize(1);
      footprint.front() = std::make_pair(oldtail, jordan.domain().lo()[0]);
    } else {  // grow
      oldtail = newtail;
      auto iter = graph.find(oldtail);
      if (iter == graph.end()) 
        throw std::runtime_error("oldtail candidates");
      cands = iter;
      if (cands->second.empty()) 
        throw std::runtime_error("cands candidates");
      outEdge = std::min_element(cands->second.cbegin(), cands->second.cend(),
                                 Selector(tol, jordan, allCrvs));
      // tmp = *outEdge;
    }
    jordan.concat(allCrvs[*outEdge]);
    newtail = jordan.endpoint();
    // delete the already-used out-edge
    removeEdge(oldtail, cands, outEdge);

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

template <int Order, typename Num, class Selector>
OPTNONE_FUNC
void PastingMap<Order, Num, Selector>::formCellClosedLoops(vector<OrientedJordanCurve<DIM, Order>> &outCont) {
  Crv jordan;
  Num oldtail;
  Num newtail;
  vector<std::pair<Num, Num>> footprint;
  auto periodCompare = [this](Num lhs, Num rhs) {
    Num dist = rhs - lhs;
    if (std::fabs(dist) < tol || std::fabs(dist + period) < tol ||
        std::fabs(dist - period) < tol)
      return 0;
    return -1;
  };

  while (!necessaryEdge.empty() || !jordan.empty()) {
    Iter outEdge;
    size_t tmp;
    auto cands = cellGraph.begin();
    if (jordan.empty()) {  // start a new loop
      auto id = *necessaryEdge.begin();
      auto startpoint = endPoints[id].first;
      auto iter = cellGraph.find(startpoint);
      cands = iter;
      outEdge = cands->second.find(id);
      // take care of the starting point
      oldtail = iter->first;
      footprint.resize(1);
      footprint.front() = std::make_pair(oldtail, jordan.domain().lo()[0]);
    } else {  // grow
      oldtail = newtail;
      auto iter = cellGraph.find(oldtail);
      if (iter == cellGraph.end()) iter = cellGraph.find(oldtail - period);
      if (iter == cellGraph.end()) iter = cellGraph.find(oldtail + period);
      if (iter == cellGraph.end()) 
        throw std::runtime_error("oldtail candidates, may enlarge Tol value.");
      cands = iter;
      if (cands->second.empty()) 
        throw std::runtime_error("cands candidates, may enlarge Tol value.");
      outEdge = std::min_element(cands->second.cbegin(), cands->second.cend(),
                                 Selector(tol, jordan, allCrvs));
      tmp = *outEdge;
    }
    jordan.concat(allCrvs[*outEdge]);
    newtail = endPoints[*outEdge].second;
    // delete the already-used out-edge
    removeEdge(oldtail, cands, outEdge);

    auto h = std::find_if(footprint.cbegin(), footprint.cend(),
                          [&](const auto &obj) {
                            return periodCompare(newtail, obj.first) == 0;
                          });
    if (h != footprint.cend()) {  // a cycle is formed
      outCont.push_back(
          jordan.extract(h->second, jordan.domain().hi()[0], tol, true));
      jordan = jordan.extract(jordan.domain().lo()[0], h->second, tol, true);
      if (!jordan.empty() &&
          jordan.getKnots().back() - jordan.getKnots().front() < tol)
        jordan = Crv();
      footprint.erase(++h, footprint.cend());
    } else {
      footprint.emplace_back(newtail, jordan.domain().hi()[0]);
    }
  }
  assert(jordan.empty());
}

template <int Order, typename Num>
bool OutEdgeSelectorByKnots<Order, Num>::operator()(const size_t &lhsId,
                                               const size_t &rhsId) const {
  return compare(lhsId, rhsId) == -1;
}
template <int Order, typename Num>
OPTNONE_FUNC
int OutEdgeSelectorByKnots<Order, Num>::compare(const size_t &lhsId,
                                           const size_t &rhsId) const {
  const auto &lhs = allCrvs[lhsId];
  const auto &rhs = allCrvs[rhsId];
  // auto d1 = normalize(lhs.midpoint() - lhs.startpoint());
  // auto d2 = normalize(rhs.midpoint() - rhs.startpoint());
  auto d1 = normalize(lhs.getComparablePoint(tol, 0) = lhs.startpoint());
  auto d2 = normalize(rhs.getComparablePoint(tol, 0) = lhs.startpoint());
  if (norm(d1 - d2) < tol) return 0;
  Num s1 = cross(d1, indir);
  Num s2 = cross(d2, indir);
  Num c1 = dot(d1, indir);
  Num c2 = dot(d2, indir);
  Num angle1 = std::atan2(s1, c1);
  if (angle1 < 0) angle1 += 2 * M_PI;
  Num angle2 = std::atan2(s2, c2);
  if (angle2 < 0) angle2 += 2 * M_PI;
  return angle1 <= angle2 ? -1 : 1;
}

//==========================================================
template class PastingMap<2>;
template class PastingMap<4>;
