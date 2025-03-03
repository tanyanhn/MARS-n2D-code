#include "PastingMap.h"

#include <algorithm>

template <int Order, class Selector>
void PastingMap<Order, Selector>::removeEdge(
    const rVec &oldtail, vector<Crv> &repo,
    typename vector<Crv>::const_iterator eit) {
  assert(repo.cbegin() <= eit && eit < repo.cend());
  repo.erase(eit);
  if (repo.empty()) graph.erase(oldtail);
}

template <int Order, class Selector>
void PastingMap<Order, Selector>::formClosedLoops(vector<Crv> &outCont) {
  Crv jordan;
  rVec oldtail, newtail;
  vector<std::pair<rVec, Real>> footprint;
  VecCompare<Real, SpaceDim> vcmp(tol);

  while (!graph.empty()) {
    typename vector<Crv>::const_iterator outEdge;
    vector<Crv> *cands;
    if (jordan.empty()) {  // start a new loop
      auto iter = graph.begin();
      cands = &(iter->second);
      outEdge = cands->cbegin();
      // take care of the starting point
      oldtail = iter->first;
      footprint.resize(1);
      footprint.front() = std::make_pair(oldtail, jordan.domain().lo()[0]);
    } else {  // grow
      oldtail = newtail;
      assert(graph.find(oldtail) != graph.end());
      cands = &(graph[oldtail]);
      outEdge = std::min_element(cands->cbegin(), cands->cend(),
                                 Selector(tol, oldtail, jordan));
    }
    jordan.concat(*outEdge);
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
bool OutEdgeSelectorByKnots<Order>::operator()(const Crv &lhs, const Crv &rhs) {
  auto getNextPoint = [](const Crv &p) {
    const auto &knots = p.getKnots();
    return p.getPolys().front()(knots[1] - knots[0]);
  };
  rVec d1 = normalize(getNextPoint(lhs) - standpoint);
  rVec d2 = normalize(getNextPoint(rhs) - standpoint);
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
