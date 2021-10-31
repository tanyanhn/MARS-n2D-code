#ifndef SPADJORBUILDER_H
#define SPADJORBUILDER_H

#include <map>
#include "Core/Curve.h"
#include "Core/VecCompare.h"

template <int Order>
class OutEdgeSelectorByKnots
{
public:
  using rVec = Vec<Real,SpaceDim>;
  using Crv = Curve<SpaceDim,Order>;

  OutEdgeSelectorByKnots(Real _tol, const rVec& tail, const Crv& existing)
    : tol(_tol), standpoint(tail)
  {
    rVec lastbutone = existing.getPolys().back()[0];
    indir = normalize(tail - lastbutone);
  }

  bool operator() (const Crv &, const Crv &);

protected:
  Real tol;
  rVec standpoint;
  rVec indir;
};

//==========================================================

template <int Order, class Selector = OutEdgeSelectorByKnots<Order>>
class PastingMap
{
public:
  using rVec = Vec<Real,SpaceDim>;
  using Crv = Curve<SpaceDim,Order>;
  template <class T> using vector = std::vector<T>;

  PastingMap(Real _tol) : tol(_tol), graph(VecCompare<Real, SpaceDim>(tol)) { }
  void formClosedLoops(vector<Crv> &outCont);
  template <class T_Crv>
  void addEdge(T_Crv&& newEdge) {
    rVec start = newEdge.startpoint();
    auto res = graph.insert(std::make_pair(start, vector<Crv> {}));
    res.first->second.push_back(std::forward<T_Crv>(newEdge));
  }

protected:
  void removeEdge(const rVec& oldtail, vector<Crv> &repo, typename vector<Crv>::const_iterator eit);
  Real tol;
  std::map<rVec, vector<Crv>, VecCompare<Real,SpaceDim>> graph;
};


#endif //SPADJORBUILDER_H
