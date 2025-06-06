#pragma once

#include "Marsn2D/Elements.h"
#include "YinSet/OrientedJordanCurve.h"
#include "YinSet/YinSet.h"

namespace Marsn2D {

template <int Order>
class approxInterfaceGraph;

class InterfaceGraph {
 public:
  template <typename T>
  using vector = std::vector<T>;
  using EdgeIndex = int;
  using SmoothnessIndicator = std::pair<EdgeIndex, EdgeIndex>;

  InterfaceGraph(vector<EdgeMark>&& edges,
                 const vector<SmoothnessIndicator>& smoothConditions);

  // auto& getFitParameter() const {
  //   return std::make_tuple(edges_, trials_, circuits_);
  // }

  template <int Order>
  friend class approxInterfaceGraph;

 protected:
  // Generate trials_ and circuits_ with edges, smoothConditions.
  static void decomposeEdgeSet(
      const vector<SmoothnessIndicator>& smoothConditions,
      const vector<EdgeMark>& edges, vector<vector<EdgeIndex>>& trials,
      vector<vector<EdgeIndex>>& circuits);

  vector<vector<EdgeIndex>> trials_;
  vector<vector<EdgeIndex>> circuits_;
  vector<EdgeMark> edges_;
  // std::vector<SmoothnessIndicator> smoothConditions_;
};

template <int Order>
class approxInterfaceGraph {
 public:
  template <typename T>
  using vector = InterfaceGraph::vector<T>;
  using OrientedJordanCurve2D = OrientedJordanCurve<DIM, Order>;
  using Trial = Curve<DIM, Order>;
  using Circuit = Curve<DIM, Order>;
  using Edge = Curve<DIM, Order>;
  using EdgeIndex = InterfaceGraph::EdgeIndex;
  using SmoothnessIndicator = InterfaceGraph::SmoothnessIndicator;

  approxInterfaceGraph(vector<EdgeMark>&& edgeMarks,
                       const vector<SmoothnessIndicator>& smoothConditions,
                       vector<vector<EdgeIndex>>&& cyclesEdgesId,
                       vector<vector<size_t>>&& YinSetId = {},
                       Real tol = distTol());

  auto approxJordanCurves() const -> vector<OrientedJordanCurve2D>;

  auto approxYinSet() const -> vector<YinSet<DIM, Order>>;

  auto accessEdges()
      -> vector<std::tuple<typename vector<Edge>::iterator,
                           typename vector<EdgeMark>::iterator, int>>;

  void updateCurve();
  auto countMarks() const -> vector<int>;
  auto countLengths() const -> vector<Real>;

  enum notAKnotBoundaryLocation {
    None = 0b00,
    Left = 0b01,
    Right = 0b10,
  };

 private:
  InterfaceGraph undirectGraph_;
  vector<Edge> edges_;
  vector<int> notaKnotBoundary_;
  // vector<Edge> reverseEdges_;
  // edges[abs(index) - 1], abs(index) starting from 1,
  // and -index imply edges[-index - 1]->reverse()
  vector<vector<EdgeIndex>> cyclesEdgesId_;
  vector<vector<size_t>> yinSetId_;
  Real tol_;
};
}  // namespace Marsn2D