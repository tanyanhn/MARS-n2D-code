#pragma once

#include "Core/FitCurve.h"
#include "Marsn2D/InterfaceGraph.h"
#include "Marsn2D/MARSReadJson.hpp"
#include "Recorder/Recorder.h"

namespace Marsn2D {

struct InterfaceGraphFactory {
  using Point = Vec<Real, 2>;
  using rVec = Point;
  using SmoothnessIndicator = InterfaceGraph::SmoothnessIndicator;
  using EdgeIndex = InterfaceGraph::EdgeIndex;
  template <typename T>
  using vector = std::vector<T>;

 public:
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
  template <int Order>
  static auto createGraph41(Real hL) -> approxInterfaceGraph<Order>;
  static auto markEllipse(Point center, rVec radius, Real hL, rVec range);
  static auto markRoseCurve(Point center, Real a, Real hL, rVec range,
                            int numPetal = 3);
  static auto markSegment(Segment<DIM> seg, Real hL);
  template <int Order>
  static auto markCurve(const Curve<DIM, Order>& crv, Real hL);
};

inline auto InterfaceGraphFactory::markRoseCurve(Point center, Real a, Real hL,
                                                 rVec range, int numPetal) {
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

inline auto InterfaceGraphFactory::markEllipse(Point center, rVec radius,
                                               Real hL, rVec range) {
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

inline auto InterfaceGraphFactory::markSegment(Segment<DIM> seg, Real hL) {
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
auto InterfaceGraphFactory::markCurve(const Curve<DIM, Order>& crv, Real hL) {
  Real knotRange = crv.getKnots().back() - crv.getKnots().front();
  auto eValue = [](auto& crv, Real t) {
    return crv(t + crv.getKnots().front());
  };
  int num = std::ceil(knotRange / hL) * 2UL;
  Real step = knotRange / num;
  vector<Point> pts;
  pts.push_back(crv.startpoint());
  for (int i = 0; i < num; i++) {
    Real localStep = step;
    Real angle = step * i + localStep;
    Point local = eValue(crv, angle);
    while (norm(pts.back() - local) > hL) {
      localStep /= 2;
      angle = step * i + localStep;
      local = eValue(crv, angle);
    }
    pts.push_back(local);
  }
  if (norm(crv.endpoint() - pts.back()) < distTol()) {
    pts.pop_back();
  }
  pts.push_back(crv.endpoint());
  return pts;
}

template <int Order>
auto InterfaceGraphFactory::createEllipse(Point center, rVec radius, Real hL,
                                          rVec range) -> YinSet<2, Order> {
  auto pts = markEllipse(center, radius, hL, range);

  auto crv = fitCurve<Order>(pts, Curve<2, Order>::periodic);
  vector<OrientedJordanCurve<2, Order>> vCrv{crv};
  YinSet<2, Order> YS(SegmentedRealizableSpadjor<Order>(vCrv), distTol());
  return YS;
}

template <int Order>
auto InterfaceGraphFactory::createDiskGraph(Point center, rVec radius, Real hL,
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
    smoothConditions.emplace_back((int)edgeMarks.size(), (int)edgeMarks.size());
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
auto InterfaceGraphFactory::createRoseGraph(Point center, Real a, Real hL,
                                            vector<Real> parts, int numPetal)
    -> approxInterfaceGraph<Order> {
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
auto InterfaceGraphFactory::createGraph41(Real hL)
    -> approxInterfaceGraph<Order> {
  vector<EdgeMark> edgeMarks(16);
  vector<SmoothnessIndicator> smoothConditions;
  vector<vector<EdgeIndex>> cyclesEdgesId(11);
  vector<vector<size_t>> YinSetId(6);
  // read tangle curve from file
  auto root = std::string(ROOT_DIR);
  std::ifstream infile(root + "/test/data1/Graph4_1.input");
  auto readPts = [](auto& infile) {
    int num;
    infile.read((char*)&num, sizeof(num));
    vector<Real> xy(2UL * num);
    infile.read((char*)xy.data(), 2UL * num * sizeof(Real));
    vector<Point> pts;
    for (int i = 0; i < num; i++) {
      pts.push_back(Point{xy[2UL * i], xy[2UL * i + 1]});
    }
    return pts;
  };

  auto pts1 = readPts(infile);
  auto pts2 = readPts(infile);
  int k[4];
  infile.read((char*)k, sizeof(k) * 4UL);
  if (norm(pts1[k[0]] - pts2[k[3]]) > distTol() ||
      norm(pts1[k[1]] - pts2[k[2]]) > distTol()) {
    throw std::runtime_error("Invalid input file");
  }
  auto rotatePts = [](const auto& pts, const auto& m) {
    vector<Point> res(pts.begin() + m, pts.end());
    res.insert(res.end(), pts.begin() + 1, pts.begin() + m + 1);
    return res;
  };
  pts1 = rotatePts(pts1, k[1]), k[0] += pts1.size() - 1 - k[1], k[1] = 0;
  if (norm(pts1[k[0]] - pts2[k[3]]) > distTol() ||
      norm(pts1[k[1]] - pts2[k[2]]) > distTol()) {
    throw std::runtime_error("Invalid input file");
  }
  auto crv1 = fitCurve<Order>(pts1, Curve<2, Order>::periodic);
  auto crv2 = fitCurve<Order>(pts2, Curve<2, Order>::notAKnot);
  vector<Curve<2, Order>> splitRes;
  crv1.split(vector<Real>{crv1.getKnots()[k[1]], crv1.getKnots()[k[0]]},
             splitRes, distTol());
  edgeMarks[0] = markCurve(splitRes[0], hL);
  edgeMarks[8] = markCurve(splitRes[1], hL);
  splitRes.clear();
  crv2.split(vector<Real>{crv2.getKnots()[0], crv2.getKnots()[k[2]],
                          crv2.getKnots()[k[3]]},
             splitRes, distTol());
  edgeMarks[12] = markCurve(splitRes[0], hL);
  edgeMarks[11] = markCurve(splitRes[1], hL);
  edgeMarks[13] = markCurve(splitRes[2], hL);
  splitRes.clear();

  smoothConditions.emplace_back(13, 12);
  smoothConditions.emplace_back(12, 14);
  smoothConditions.emplace_back(1, 9);
  smoothConditions.emplace_back(9, 1);
  cyclesEdgesId[3] = (vector<EdgeIndex>{1, 9});
  cyclesEdgesId[6] = (vector<EdgeIndex>{13, 12, 14});
  YinSetId[3] = (vector<size_t>{3, 6});

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
    const auto disk = InterfaceGraphFactory::createDiskGraph<Order>(
        center, radius, hL / 2, parts);
    vecDisk.emplace_back(std::move(disk));
    vecN.emplace_back(N);
    vecBox.emplace_back(0, N - 1);
    // vecH.emplace_back(h);
    vecHL.emplace_back(hL);
    vecDt.emplace_back(dt);
    N = N * 2;
  }

  // accurate solution
  Real hL = vecHL.back() * rTiny * 10;
  auto exactDisk = InterfaceGraphFactory::createDiskGraph<Order>(center, radius,
                                                                 hL / 2, parts);

  return make_tuple(vecDisk, exactDisk, radius, exactArea, exactLength, vecBox,
                    vecN, aimOrder, vecHL, rTiny, nGrid, curvConfig, plotConfig,
                    printDetail, t0, vecDt, te, timeIntegrator, velocityPtr);
}

}  // namespace Marsn2D