#ifndef SPADJORBUILDER_H
#define SPADJORBUILDER_H

#include <map>
#include <set>

#include "Core/Curve.h"
#include "Core/VecCompare.h"
#include "YinSet/OrientedJordanCurve.h"

template <int Order, typename Num = Real>
class OutEdgeSelectorByKnots {
 public:
  using rVec = Vec<Num, SpaceDim>;
  using Crv = Curve<SpaceDim, Order>;

  OutEdgeSelectorByKnots(Num _tol, const Crv& existing,
                         const std::vector<Crv>& allCrvs_)
      : tol(_tol),
        // indir(existing.midpoint() - existing.endpoint()),
        indir(normalize(existing.getComparablePoint(_tol, 1) -
                        existing.endpoint())),
        allCrvs(allCrvs_) {}

  int compare(const size_t&, const size_t&) const;
  bool operator()(const size_t&, const size_t&) const;

 protected:
  Num tol;
  // rVec standpoint;
  rVec indir;
  const std::vector<Crv>& allCrvs;
};

//==========================================================

template <int Order, typename Num = Real,
          class Selector = OutEdgeSelectorByKnots<Order, Num>>
class PastingMap {
 public:
  using rVec = Vec<Num, SpaceDim>;
  using Crv = Curve<SpaceDim, Order>;
  template <class T>
  using vector = std::vector<T>;

  PastingMap(Num _tol)
      : tol(_tol),
        graph(VecCompare<Num, SpaceDim>(tol)),
        cellGraph(VecCompare<Num, SpaceDim - 1>(tol)) {}
  void formClosedLoops(vector<OrientedJordanCurve<DIM, Order>>& outCont);
  void formCellClosedLoops(vector<OrientedJordanCurve<DIM, Order>>& outCont);
  template <class T_Crv>
  void addEdge(T_Crv&& newEdge, bool necessary = true);
  template <class T_Crv>
  void addCellEdge(T_Crv&& newEdge, Num start, Num end,
                   bool necessary = true);
  void setPeriod(Num p) { period = p; }
  void addVertex(const rVec& newVertex) { graph.insert({newVertex, {}}); }

 protected:
  using Iter = std::set<size_t>::const_iterator;
  void removeEdge(const rVec& oldtail,  auto& repo, Iter eit);
  Num tol;
  std::vector<Crv> allCrvs;
  std::map<rVec, std::set<size_t>, VecCompare<Num, SpaceDim>> graph;
  std::map<Num, std::set<size_t>, VecCompare<Num, 1>> cellGraph;
  std::set<size_t> necessaryEdge;  // for the necessary edge form loop
  bool useCellGraph = false;
  std::unordered_map<size_t, std::pair<Num, Num>> endPoints;
  Num period;
};

template <int Order, typename Num, class Selector>
template <class T_Crv>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
void PastingMap<Order, Num, Selector>::addEdge(T_Crv&& newEdge, bool necessary) {
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


template <int Order, typename Num, class Selector>
template <class T_Crv>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
void PastingMap<Order, Num, Selector>::addCellEdge(T_Crv&& newEdge, Num start, Num end, bool necessary) {
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
