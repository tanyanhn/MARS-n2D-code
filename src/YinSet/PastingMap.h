#ifndef SPADJORBUILDER_H
#define SPADJORBUILDER_H

#include <map>
#include <set>

#include "Core/Curve.h"
#include "Core/VecCompare.h"
#include "YinSet/OrientedJordanCurve.h"

template <int Order>
class OutEdgeSelectorByKnots {
 public:
  using rVec = Vec<Real, SpaceDim>;
  using Crv = Curve<SpaceDim, Order>;

  OutEdgeSelectorByKnots(Real _tol, const Crv& existing,
                         const std::vector<Crv>& allCrvs_)
      : tol(_tol),
        indir(existing.midpoint() - existing.endpoint()),
        // indir(existing.getComparableDirection(_tol, 1)),
        allCrvs(allCrvs_) {}

  int compare(const size_t&, const size_t&) const;
  bool operator()(const size_t&, const size_t&) const;

 protected:
  Real tol;
  // rVec standpoint;
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

  PastingMap(Real _tol)
      : tol(_tol),
        graph(VecCompare<Real, SpaceDim>(tol)),
        cellGraph(VecCompare<Real, SpaceDim - 1>(tol)) {}
  void formClosedLoops(vector<OrientedJordanCurve<DIM, Order>>& outCont);
  void formCellClosedLoops(vector<OrientedJordanCurve<DIM, Order>>& outCont);
  template <class T_Crv>
  void addEdge(T_Crv&& newEdge, bool necessary = true);
  template <class T_Crv>
  void addCellEdge(T_Crv&& newEdge, Real start, Real end,
                   bool necessary = true);
  void setPeriod(Real p) { period = p; }
  void addVertex(const rVec& newVertex) { graph.insert({newVertex, {}}); }

 protected:
  using Iter = std::set<size_t>::const_iterator;
  void removeEdge(const rVec& oldtail,  auto& repo, Iter eit);
  Real tol;
  std::vector<Crv> allCrvs;
  std::map<rVec, std::set<size_t>, VecCompare<Real, SpaceDim>> graph;
  std::map<Real, std::set<size_t>, VecCompare<Real, 1>> cellGraph;
  std::set<size_t> necessaryEdge;  // for the necessary edge form loop
  bool useCellGraph = false;
  std::unordered_map<size_t, std::pair<Real, Real>> endPoints;
  Real period;
};

template <int Order, class Selector>
template <class T_Crv>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
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
  // remove overlap reverse edge.
  auto reverseCrv = newEdge.reverse();
  iter = graph.find(reverseCrv.startpoint());
  if (iter != graph.end()) {
    for (auto id : iter->second) {
      if (reverseCrv.equal(allCrvs[id], tol)) {
        necessaryEdge.erase(id);
        iter->second.erase(id);
        if (iter->second.empty()) graph.erase(iter);
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


template <int Order, class Selector>
template <class T_Crv>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
void PastingMap<Order, Selector>::addCellEdge(T_Crv&& newEdge, Real start, Real end, bool necessary) {
  useCellGraph = true;
  auto iter = cellGraph.find(start);
  if (iter != cellGraph.end()) {
    for (auto id : iter->second) {
      if (newEdge.equal(allCrvs[id], tol)) {
        if (necessary) {
          necessaryEdge.insert(id);
        }
        return;
      }
    }
  }
  // remove overlap reverse edge.
  auto reverseCrv = newEdge.reverse();
  iter = cellGraph.find(end);
  if (iter != cellGraph.end()) {
    for (auto id : iter->second) {
      if (reverseCrv.equal(allCrvs[id], tol)) {
        necessaryEdge.erase(id);
        iter->second.erase(id);
        if (iter->second.empty()) cellGraph.erase(iter);
        return;
      }
    }
  }
  allCrvs.push_back(std::forward<T_Crv>(newEdge));
  size_t id = allCrvs.size() - 1;
  auto res = cellGraph.insert(std::make_pair(start, std::set<size_t>{}));
  res.first->second.insert(id);
  endPoints[id] = std::make_pair(start, end);
  if (necessary) necessaryEdge.insert(id);
}

#endif  // SPADJORBUILDER_H
