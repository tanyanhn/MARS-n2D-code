#pragma once

#include <memory>

#include "Core/Curve.h"
#include "Core/Vec.h"

namespace MARSn2D {
using Vertex = Vec<Real, DIM>;

using Junction = Vec<Real, DIM>;

using Kink = Vec<Real, DIM>;

using Basepoint = Vec<Real, DIM>;

template <int Order>
using Edge = Curve<DIM, Order>;

using EdgeMark = std::vector<Vertex>;

template <int Order>
using EdgePtr = std::shared_ptr<Edge<Order>>;

using EdgeMarkPtr = std::shared_ptr<EdgeMark>;

template <int Order>
using Trial = Curve<DIM, Order>;

template <int Order>
using Circuit = Curve<DIM, Order>;

template <int Order>
using Cycle = Curve<DIM, Order>;

// template <int Order>
// using SmoothnessIndicator = std::pair<EdgePtr<Order>, EdgePtr<Order>>;

}  // namespace MARSn2D
