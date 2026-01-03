#pragma once
#ifndef HLGENERATORFACTORY_H
#define HLGENERATORFACTORY_H

#include <functional>

#include <Core/Config.h>

namespace Marsn2D {

struct HLGeneratorFactory {
  using Prod =
      std::function<void(const std::vector<Real>& curv, std::vector<Real>& hL,
                         Real& lowerScale, Real& upperScale, Real hL_)>;
  static auto create(Real rCMin, Real rhoMin, Real rhoMax, Real rTiny,
              const std::function<Real(Real)>& sigma) -> Prod;
};

inline auto HLGeneratorFactory::create(Real rCMin, Real rhoMin, Real rhoMax,
                                       Real rTiny,
                                       const std::function<Real(Real)>& sigma)
    -> Prod {
  return [rCMin, rhoMin, rhoMax, rTiny, sigma](
             const std::vector<Real>& curv, std::vector<Real>& hL,
             Real& lowerScale, Real& upperScale, Real hL_) {
    lowerScale = rTiny;
    upperScale = 1 - 2 * rTiny;

    hL.resize(curv.size());
    Real curvMax = curv.front();
    Real curvMin = curv.front();
    for (const auto& c : curv) {
      curvMax = std::max(curvMax, c);
      curvMin = std::min(curvMin, c);
    }
    curvMax = std::min(curvMax, rhoMax);
    curvMin = std::max(curvMin, rhoMin);
    Real rMin = std::max(rCMin, curvMin / curvMax);

    for (size_t i = 0; i < curv.size(); i++) {
      if (curv[i] <= curvMin)
        hL[i] = hL_;
      else if (curv[i] >= curvMax)
        hL[i] = hL_ * rMin;
      else {
        hL[i] = rMin * hL_ + (1 - rMin) * hL_ *
                                 sigma((1 / curv[i] - 1 / curvMax) /
                                       (1 / curvMin - 1 / curvMax));
      }
    }
  };
}
}  // namespace Marsn2D
#endif // !HLGENERATORFACTORY_H
