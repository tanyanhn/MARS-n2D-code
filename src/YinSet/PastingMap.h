#ifndef SPADJORBUILDER_H
#define SPADJORBUILDER_H

#include <map>
#include <set>

#include "Core/Curve.h"
#include "Core/VecCompare.h"

template <int Order>
class OutEdgeSelectorByKnots {
 public:
  using rVec = Vec<Real, SpaceDim>;
  using Crv = Curve<SpaceDim, Order>;

  OutEdgeSelectorByKnots(Real _tol, const Crv& existing,
                         const std::vector<Crv>& allCrvs_)
      : tol(_tol), standpoint(existing.endpoint()), allCrvs(allCrvs_) {
    rVec lastButOne = existing.getComparablePoint(_tol, 1);
    indir = normalize(lastButOne - standpoint);
  }

  int compare(const size_t&, const size_t&) const;
  bool operator()(const size_t&, const size_t&) const;

 protected:
  Real tol;
  rVec standpoint;
  rVec indir;
  const std::vector<Crv>& allCrvs;
};

//==========================================================

template <int Order, class Selector = OutEdgeSelectorByKnots<Order>>
class PastingMap {
 public:
  using rVec = Vec<Real, SpaceDim>;
  using Crv = Curve<SpaceDim, Order>;
  template <class T>
  using vector = std::vector<T>;

  PastingMap(Real _tol) : tol(_tol), graph(VecCompare<Real, SpaceDim>(tol)) {}
  void formClosedLoops(vector<Crv>& outCont);
  template <class T_Crv>
  void addEdge(T_Crv&& newEdge, bool necessary = true);

 protected:
  using Iter = std::set<size_t>::const_iterator;
  void removeEdge(const rVec& oldtail, std::set<size_t>& repo, Iter eit);
  Real tol;
  std::vector<Crv> allCrvs;
  std::map<rVec, std::set<size_t>, VecCompare<Real, SpaceDim>> graph;
  std::set<size_t> necessaryEdge;  // for the necessary edge form loop
};

template <int Order, class Selector>
template <class T_Crv>
void PastingMap<Order, Selector>::addEdge(T_Crv&& newEdge, bool necessary) {
  rVec start = newEdge.startpoint();
  auto iter = graph.find(start);
  if (iter != graph.end()) {
    for (auto id : iter->second) {
      if (newEdge.equal(allCrvs[id], tol)) {
        if (necessary) {
          necessaryEdge.insert(id);
        }
        return;
      }
    }
  }
  allCrvs.push_back(std::forward<T_Crv>(newEdge));
  size_t id = allCrvs.size() - 1;
  auto res = graph.insert(std::make_pair(start, std::set<size_t>{}));
  res.first->second.insert(id);
  if (necessary) necessaryEdge.insert(id);
}

#endif  // SPADJORBUILDER_H
