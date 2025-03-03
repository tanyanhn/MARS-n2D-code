#pragma once

#include "Mars/Elements.h"
#include "Mars/OrientedJordanCurve.h"

namespace MARSn2D {

template <int Order>
class approxInterfaceGraph;

class InterfaceGraph {
 public:
  using EdgeIndex = int;
  using SmoothnessIndicator = std::pair<EdgeIndex, EdgeIndex>;

  InterfaceGraph(std::vector<EdgeMark>&& edges,
                 std::vector<SmoothnessIndicator>&& smoothConditions);

  // auto& getFitParameter() const {
  //   return std::make_tuple(edges_, trials_, circuits_);
  // }

  template <int Order>
  friend class approxInterfaceGraph;

 protected:
  // Generate trials_ and circuits_ with edges, smoothConditions.
  static void decomposeEdgeSet(const auto& edges, const auto& smoothConditions,
                               auto& trials, auto& circuits);

  std::vector<std::vector<EdgeIndex>> trials_;
  std::vector<std::vector<EdgeIndex>> circuits_;
  std::vector<EdgeMark> edges_;
  // std::vector<SmoothnessIndicator> smoothConditions_;
};

template <int Order>
class approxInterfaceGraph {
 public:
  using OrientedJordanCurve = OrientedJordanCurve<Order>;
  using Trial = Trial<Order>;
  using Circuit = Circuit<Order>;
  using Edge = Edge<Order>;
  using EdgeIndex = InterfaceGraph::EdgeIndex;
  using SmoothnessIndicator = InterfaceGraph::SmoothnessIndicator;

  approxInterfaceGraph(std::vector<EdgeMark>&& edges,
                       std::vector<SmoothnessIndicator>&& smoothConditions,
                       // edges[abs(index) - 1], abs(index) starting from 1,
                       // and -index imply edges[-index - 1]->reverse()
                       std::vector<std::vector<EdgeIndex>>&& cyclesEdgesId);

  auto& approxJordanCurves() const;

 private:
  InterfaceGraph undirectGraph;
  std::vector<Edge> edges_;
  std::vector<std::vector<EdgeIndex>>& cyclesEdgesId_;
};
}  // namespace MARSn2D