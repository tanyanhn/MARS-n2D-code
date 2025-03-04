#pragma once

#include <unordered_set>

#include "Marsn2D/Elements.h"

namespace MARSn2D {

template <int Order>
class OrientedJordanCurve {
 public:
  using Edge = Edge<Order>;

 private:
  std::vector<Edge> cycle_;
  std::unordered_set<size_t> kinks_;
  std::unordered_set<size_t> junctions_;
  // basepoint always on index 0
  // size_t basepoint = 0;
};
}  // namespace MARSn2D