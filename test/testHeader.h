#pragma once
// Only include by test file.
#include "Core/Config.h"
#include "Core/Vec.h"
#include "Core/dirConfig.h"
#include "Marsn2D/InterfaceGraph.h"
#include "Marsn2D/InterfaceGraphFactory.h"
#include "Marsn2D/MARSReadJson.hpp"
#include "catch_amalgamated.hpp"

using namespace std;
using namespace Catch;
using namespace Marsn2D;

using Point = Vec<Real, 2>;
using rVec = Point;

const std::string rootDir(ROOT_DIR);

struct Generator {
  static VecCompare<Real, 2> vCmp;
  template <int l, int r>
  [[nodiscard]] static auto randomCreateReal() -> Real;

  template <int Order>
  static auto createEllipse(Point center, rVec radius, Real hL,
                            rVec range = {0, 2 * M_PI}) -> YinSet<2, Order>;

  template <int Order>
  static auto createDiskGraph(Point center, rVec radius, Real hL,
                              vector<Real> parts)
      -> approxInterfaceGraph<Order>;
  template <int Order>
  static auto createRoseGraph(Point center, Real a, Real hL, vector<Real> parts,
                              int numPetal = 3) -> approxInterfaceGraph<Order>;
  static auto markEllipse(Point center, rVec radius, Real hL, rVec range);
  static auto markRoseCurve(Point center, Real a, Real hL, rVec range,
                            int numPetal = 3);
  static auto markSegment(Segment<DIM> seg, Real hL);
};

const Generator generator;
VecCompare<Real, 2> Generator::vCmp(distTol());

template <int l, int r>
auto Generator::randomCreateReal() -> Real {
  STATIC_CHECK(l < r);
  return GENERATE(take(1, random(l, r)));
}

inline auto Generator::markRoseCurve(Point center, Real a, Real hL, rVec range,
                                     int numPetal) {
  return InterfaceGraphFactory::markRoseCurve(center, a, hL, range, numPetal);
}

inline auto Generator::markEllipse(Point center, rVec radius, Real hL,
                                   rVec range) {
  return InterfaceGraphFactory::markEllipse(center, radius, hL, range);
}

inline auto Generator::markSegment(Segment<DIM> seg, Real hL) {
  return InterfaceGraphFactory::markSegment(seg, hL);
}

template <int Order>
auto Generator::createEllipse(Point center, rVec radius, Real hL, rVec range)
    -> YinSet<2, Order> {
  return InterfaceGraphFactory::createEllipse<Order>(center, radius, hL, range);
}

template <int Order>
auto Generator::createDiskGraph(Point center, rVec radius, Real hL,
                                vector<Real> parts)
    -> approxInterfaceGraph<Order> {
  return InterfaceGraphFactory::createDiskGraph<Order>(center, radius, hL,
                                                       parts);
}

template <int Order>
auto Generator::createRoseGraph(Point center, Real a, Real hL,
                                vector<Real> parts, int numPetal)
    -> approxInterfaceGraph<Order> {
  return InterfaceGraphFactory::createRoseGraph<Order>(center, a, hL, parts,
                                                       numPetal);
}
