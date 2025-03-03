#pragma once

#include "Mars/InterfaceGraph.h"
#include "Mars/Yinset.h"

namespace MARSn2D {
template <int Order>
class InterfaceTrack {
 public:
 struct curvBaseStrategy {
   Real rhoMin;
   Real rhoMax;
   Real rCMin;
   bool curvH = false;
 };


  std::vector<Yinset<Order>> operator()(
      const std::vector<Yinset<Order>>& input, Real t0, Real te, auto vel,
      Real rTiny, Real hL, int OrderInTime,
      curvBaseStrategy curv);
      
  auto ARMSStep(Real ts0, Real tse, auto vel, Real rTiny, Real hL,
                int OrderInTime);
};
}  // namespace MARSn2D