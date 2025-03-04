#pragma once

#include <vector>

#include "Marsn2D/Elements.h"
#include "Marsn2D/OrientedJordanCurve.h"

namespace MARSn2D {
template <int Order>
class Spadjor {
 public:
  using OrientedJordanCurve = OrientedJordanCurve<Order>;

 protected:
  std::vector<OrientedJordanCurve> orientedJordanCurves_;
};

template <int Order>
class Yinset : public Spadjor<Order> {
 public:
  /// A node in the Hasse diagram.
  struct Node {
    int depth;  // even number for positive orientation, odd number of negative
                // orientation
    int parent;
    std::vector<int> children;
  };

 protected:
  void buildHasse(Real tol);

  /// The Hasse diagram. The last node is the root of the forest.
  std::vector<Node> diagram_;
};
}  // namespace MARSn2D