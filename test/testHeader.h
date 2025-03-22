#pragma once
// Only include by test file.
#include "Core/Config.h"
#include "Core/Vec.h"
#include "Core/dirConfig.h"
#include "Marsn2D/InterfaceGraph.h"
#include "Marsn2D/MARSReadJson.hpp"
#include "YinSet/YinSet.h"
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
  vector<Point> pts;
  VecCompare<Real, DIM> vCmp(distTol());
  // |x'(\theta)| \leq 4 a
  // \sqrt{(x(\theta1) - x(\theta2))^2 + (y(\theta1) - y(\theta2))^2} \leq 4
  // \sqrt{2} a \Delta \theta.
  Real stepAngle = hL / a / 4 / sqrt(2);
  int step = std::ceil((range[1] - range[0]) / stepAngle);
  stepAngle = (range[1] - range[0]) / step;
  auto eValue = [a, center, range, numPetal](Real angle) {
    angle += range[0];
    Real radius = a * sin(numPetal * angle);
    return Point{center + rVec{radius * cos(angle), radius * sin(angle)}};
  };
  pts.reserve(step + 1);
  pts.push_back(eValue(0));
  for (int i = 0; i < step; i++) {
    Real localAngle = stepAngle;
    Real angle = stepAngle * i + localAngle;
    Point local = eValue(angle);
    while (norm(pts.back() - local) > hL) {
      localAngle /= 2;
      angle = stepAngle * i + localAngle;
      local = eValue(angle);
    }

    while (angle < stepAngle * (i + 1) - distTol()) {
      pts.push_back(eValue(angle));
      angle += localAngle;
    }
    pts.push_back(eValue(stepAngle * (i + 1)));
  }
  auto p1 = eValue(range[1] - range[0]);
  if (vCmp(pts.back(), p1) == 0) pts.pop_back();

  pts.push_back(p1);
  return pts;
}

inline auto Generator::markEllipse(Point center, rVec radius, Real hL,
                                   rVec range) {
  vector<Point> pts;
  VecCompare<Real, DIM> vCmp(distTol());
  Real stepAngle = hL / std::sqrt(radius[0] * radius[1]);
  if (range[1] < range[0]) range[1] += 2 * M_PI;
  int step = std::ceil((range[1] - range[0]) / stepAngle);
  stepAngle = (range[1] - range[0]) / step;
  auto eValue = [radius, center, range](Real angle) {
    angle += range[0];
    return Point{center + rVec{radius[0] * cos(angle), radius[1] * sin(angle)}};
  };
  pts.reserve(step + 1);
  pts.push_back(eValue(0));
  for (int i = 0; i < step; i++) {
    Real localAngle = stepAngle;
    Real angle = stepAngle * i + localAngle;
    Point local = eValue(angle);
    while (norm(pts.back() - local) > hL) {
      localAngle /= 2;
      angle = stepAngle * i + localAngle;
      local = eValue(angle);
    }

    while (angle < stepAngle * (i + 1) - distTol()) {
      pts.push_back(eValue(angle));
      angle += localAngle;
    }
    pts.push_back(eValue(stepAngle * (i + 1)));
  }
  auto p1 = eValue(range[1] - range[0]);
  if (vCmp(pts.back(), p1) == 0) pts.pop_back();

  pts.push_back(p1);
  return pts;
}

inline auto Generator::markSegment(Segment<DIM> seg, Real hL) {
  vector<Point> pts;
  VecCompare<Real, DIM> vCmp(distTol());
  auto p0 = seg[0];
  auto p1 = seg[1];
  Real length = norm(p1 - p0);
  int num = std::ceil(length / hL);
  rVec piece = (p1 - p0) / num;
  Point loc = p0;
  for (int i = 0; i < num; i++) {
    pts.push_back(loc);
    loc = p0 + piece * (i + 1);
  }

  if (vCmp(pts.back(), p1) == 0) pts.pop_back();
  pts.push_back(p1);

  return pts;
}

template <int Order>
auto Generator::createEllipse(Point center, rVec radius, Real hL, rVec range)
    -> YinSet<2, Order> {
  auto pts = markEllipse(center, radius, hL, range);

  auto crv = fitCurve<Order>(pts, Curve<2, Order>::periodic);
  vector<OrientedJordanCurve<2, Order>> vCrv{crv};
  YinSet<2, Order> YS(SegmentedRealizableSpadjor<Order>(vCrv), distTol());
  return YS;
}

template <int Order>
auto Generator::createDiskGraph(Point center, rVec radius, Real hL,
                                vector<Real> parts)
    -> approxInterfaceGraph<Order> {
  using SmoothnessIndicator = InterfaceGraph::SmoothnessIndicator;
  using EdgeIndex = InterfaceGraph::EdgeIndex;
  using std::vector;

  vector<EdgeMark> edgeMarks;
  vector<SmoothnessIndicator> smoothConditions;
  vector<vector<EdgeIndex>> cyclesEdgesId(parts.size() + 1);
  vector<vector<size_t>> YinSetId(parts.size() + 1);

  if (parts.empty()) {
    cyclesEdgesId.resize(2);
    YinSetId.resize(2);
    EdgeMark pts = markEllipse(center, radius, hL, {0, 2 * M_PI});
    edgeMarks.push_back(std::move(pts));
    smoothConditions.emplace_back((int)edgeMarks.size(),
                                  (int)edgeMarks.size());
    cyclesEdgesId[0].push_back(1);
    cyclesEdgesId[1].push_back(-1);
    YinSetId[0].push_back(0);
    YinSetId[1].push_back(1);
  } else {
    auto eValue = [center, radius](Real angle) {
      return Point{center +
                   rVec{radius[0] * cos(angle), radius[1] * sin(angle)}};
    };

    vector<Segment<DIM>> separateSegment;
    for (size_t i = 0; i < parts.size(); ++i) {
      Segment seg(eValue(parts[i]), center);
      int parallelSegId = -1;
      for (size_t j = 0; j < separateSegment.size(); ++j) {
        if (seg.parallel(separateSegment[j], distTol())) {
          parallelSegId = j;
          break;
        }
      }
      separateSegment.push_back(seg);
      auto pts = markSegment(seg, hL);
      int sign = -2;
      if (parallelSegId != -1) {
        std::reverse(pts.begin(), pts.end());
        smoothConditions.emplace_back(parallelSegId + 1,
                                      (int)edgeMarks.size() + 1);
        sign = -1;
      } else {
        sign = 1;
      }
      edgeMarks.push_back(std::move(pts));
      auto& vec = cyclesEdgesId[(i - 1 + parts.size()) % parts.size()];
      vec.push_back(sign * (int)edgeMarks.size());
      if (vec.size() == 2) std::swap(vec[0], vec[1]);
      cyclesEdgesId[i].push_back(-sign * (int)edgeMarks.size());
    }

    if (norm(parts.back() - parts.front()) > distTol())
      parts.push_back(parts.front());

    for (size_t i = 0; i < parts.size() - 1; ++i) {
      EdgeMark pts = markEllipse(center, radius, hL, {parts[i], parts[i + 1]});

      edgeMarks.push_back(std::move(pts));
      smoothConditions.emplace_back(
          (int)edgeMarks.size() - 1 + (i == 0 ? parts.size() - 1 : 0),
          (int)edgeMarks.size());
      cyclesEdgesId[i].push_back((int)edgeMarks.size());
      cyclesEdgesId[parts.size() - 1].push_back(-(int)edgeMarks.size());
      YinSetId[i].push_back(i);
    }
    std::reverse(cyclesEdgesId[parts.size() - 1].begin(),
                 cyclesEdgesId[parts.size() - 1].end());
    YinSetId[parts.size() - 1].push_back(parts.size() - 1);
    // smoothConditions.clear();
  }

  return approxInterfaceGraph<Order>(std::move(edgeMarks), smoothConditions,
                                     std::move(cyclesEdgesId),
                                     std::move(YinSetId), distTol());
}

template <int Order>
auto Generator::createRoseGraph(Point center, Real a, Real hL,
                                vector<Real> parts, int numPetal)
    -> approxInterfaceGraph<Order> {
  using SmoothnessIndicator = InterfaceGraph::SmoothnessIndicator;
  using EdgeIndex = InterfaceGraph::EdgeIndex;
  using std::vector;

  vector<EdgeMark> edgeMarks;
  vector<SmoothnessIndicator> smoothConditions;
  vector<vector<EdgeIndex>> cyclesEdgesId(parts.size());
  vector<vector<size_t>> YinSetId(parts.size());

  // r = a sin(n \theta)
  // x = r cos(\theta), y = r sin(\theta)
  // auto eValue = [center, a, numPetal](Real angle) {
  //   Real r = a * sin(numPetal * angle);
  //   return Point{center + rVec{r * cos(angle), r * sin(angle)}};
  // };
  Real halfAngle = M_PI / numPetal / 2;
  for (size_t i = 0; i < parts.size() - 1; ++i) {
    for (size_t j = 0; j < 2; ++j) {  // each petal is slitted into two parts
      EdgeMark pts;
      if (j == 0) {
        pts = markRoseCurve(center, a, hL, {parts[i], parts[i] + halfAngle});
        YinSetId[i].push_back(i);
      } else {
        pts = markRoseCurve(center, a, hL,
                            {parts[i + 1] - halfAngle, parts[i + 1]});
      }
      edgeMarks.push_back(std::move(pts));
      cyclesEdgesId[i].push_back((int)edgeMarks.size());
      cyclesEdgesId[parts.size() - 1].push_back(-(int)edgeMarks.size());
      if (parts.size() == numPetal + 1 && i * 2 + j == parts.size() * 2 - 3) {
        smoothConditions.emplace_back((int)edgeMarks.size(), 1);
      } else if (i * 2 + j != parts.size() * 2 - 3) {
        smoothConditions.emplace_back((int)edgeMarks.size(),
                                      (int)edgeMarks.size() + 1);
      }
    }
  }
  std::reverse(cyclesEdgesId[parts.size() - 1].begin(),
               cyclesEdgesId[parts.size() - 1].end());
  YinSetId[parts.size() - 1].push_back(parts.size() - 1);

  return approxInterfaceGraph<Order>(std::move(edgeMarks), smoothConditions,
                                     std::move(cyclesEdgesId),
                                     std::move(YinSetId), distTol());
}

template <int Order>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
auto diskTEST(const std::string& jsonFile) {
  using namespace std;
  auto params = read_config(jsonFile);
  distTol(params.tolerance.distTol);
  newtonTol(params.tolerance.newtonTol);
  newtonMaxIter(params.tolerance.newtonMaxIter);
  // time
  const Real te = params.time.te;
  const Real t0 = params.time.t0;
  const Real Cr = params.time.Cr;
  const Real uM = params.time.uM;
  // grid
  const int N0 = params.grid.N0;
  const rVec lo = params.grid.lo;
  const rVec hi = params.grid.hi;
  const int aimOrder = params.grid.aimOrder;
  const Real hLCoefficient = params.grid.hLCoefficient;
  const Real rTiny = params.grid.rTiny;
  const int nGrid = params.grid.nGrid;
  std::shared_ptr<TimeIntegrator<DIM, VectorFunction>> timeIntegrator;
  if (aimOrder == 4) {
    timeIntegrator = make_shared<ERK<DIM, RK::ClassicRK4, VectorFunction>>();
  } else if (aimOrder == 6) {
    timeIntegrator = make_shared<ERK<DIM, RK::Verner6, VectorFunction>>();
  } else if (aimOrder == 8) {
    timeIntegrator =
        make_shared<ERK<DIM, RK::PrinceDormand8, VectorFunction>>();
  }
  const Real T = te;
  const auto velocityPtr = params.field.velocity;

  // curvature-based adaption
  const bool curvUsed = params.curvature.curvUsed;
  const Real rCMin = params.curvature.rCMin;
  const Real rhoMin = params.curvature.rhoMin;
  const Real rhoMax = params.curvature.rhoMax;
  const typename Marsn2D::MARSn2D<
      Order, VectorFunction>::CurvatureAdaptionConfig curvConfig{
      .used = curvUsed,
      .rCMin = rCMin,
      .rhoCRange = Interval<1>{rhoMin, rhoMax},
      .sigma = [](Real x) { return x; },
  };

  // plot
  const int plotOutput = params.plotConfig.plotOutput;
  const bool plotInner = params.plotConfig.plotInner;
  const int plotStride = params.plotConfig.plotStride;
  const bool printDetail = params.plotConfig.printDetail;
  const struct Marsn2D::MARSn2D<Order, VectorFunction>::PlotConfig plotConfig {
    .output = plotOutput, .plotInner = plotInner, .opStride = plotStride,
    .fName = to_string(Order) + "Circle" + "_grid" + to_string(N0) + "_T" +
             to_string((int)T),
    .range = {lo, hi},
  };

  // domain disk
  const Point center = params.domain.center;
  const rVec radius = params.domain.radius;
  vector<Real> parts = params.domain.parts;
  for (auto& part : parts) part *= 2 * M_PI;
  const Real exactArea = M_PI * radius[0] * radius[1] / parts.size();
  const Real exactLength = 2 * M_PI * radius[0] / parts.size() + 2 * radius[0];

  // for multi grid.
  vector<Marsn2D::approxInterfaceGraph<Order>> vecDisk;
  vector<int> vecN;
  vector<Box<DIM>> vecBox;
  // vector<rVec> vecH;
  vector<Real> vecHL;
  vector<Real> vecDt;
  Real N = N0;
  for (uint i = 0; i < nGrid; ++i) {
    rVec h = (hi - lo) / N;
    Real hL = hLCoefficient * std::pow(h[0] * h[1], 0.5 * aimOrder / Order);
    Real dt = std::min(h[0], h[1]) / uM * Cr;
    const auto disk =
        Generator::createDiskGraph<Order>(center, radius, hL / 2, parts);
    vecDisk.emplace_back(std::move(disk));
    vecN.emplace_back(N);
    vecBox.emplace_back(0, N - 1);
    // vecH.emplace_back(h);
    vecHL.emplace_back(hL);
    vecDt.emplace_back(dt);
    N = N * 2;
  }

  // accurate solution
  Real hL = vecHL.back() / 10;
  auto exactDisk =
      Generator::createDiskGraph<Order>(center, radius, hL / 2, parts);

  return make_tuple(vecDisk, exactDisk, radius, exactArea, exactLength, vecBox,
                    vecN, aimOrder, vecHL, rTiny, nGrid, curvConfig, plotConfig,
                    printDetail, t0, vecDt, te, timeIntegrator, velocityPtr);
}