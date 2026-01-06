#pragma once

#include "Core/FitCurve.h"
#include "Marsn2D/InterfaceGraph.h"
#include "Marsn2D/MARSReadJson.hpp"
#include "Recorder/Recorder.h"
#include <tuple>
#include <fstream>
#include <algorithm>

namespace Marsn2D {

struct InterfaceGraphFactory {
  using Point = Vec<Real, 2>;
  using rVec = Point;
  using SmoothnessIndicator = InterfaceGraph::SmoothnessIndicator;
  using EdgeIndex = InterfaceGraph::EdgeIndex;
  template <typename T>
  using vector = std::vector<T>;

  template <int Order>
  using MARScurvConfig =
      typename Marsn2D::MARSn2D<Order, VectorFunction>::CurvatureAdaptionConfig;

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
  static auto markCurve(const Curve<DIM, Order>& crv, Real hL, Real lo = 1,
                        Real hi = -1);
  template <int Order>
  static auto markCurve(const Curve<DIM, Order>& crv, Real hL, Real t0,
                        Real t1, const MARScurvConfig<Order-2>& curvConfig);
  struct SplineCurvePP {
    uint32_t order{};
    std::vector<Real> breaks;
    std::vector<Real> coefsX;
    std::vector<Real> coefsY;
  };
  struct SplineEdgeRef {
    int v1{}, v2{}, idx{};
    int curveId{-1};
    int dir{0};   // +1 v1->v2, -1 reverse
    Real t0{0}, t1{0};
  };
  struct BoundaryEdge {
    int edgeId{};
    int dir{};  // +1 v1->v2, -1 reverse
  };
  struct BoundaryLoop {
    std::array<int, 3> color{};
    vector<BoundaryEdge> edges;
  };
  struct SmoothPair {
    int eA{}, fA{}, eB{}, fB{};
  };
  template <int Order>
  static auto inputFromSVG(const std::string& binPath, Real hL,
                           const MARScurvConfig<Order>& curvConfig = {})
      -> approxInterfaceGraph<Order>;
  template <int Order>
  static auto createEmpty()
      -> approxInterfaceGraph<Order>;
  template <int Order>
  static auto dealRaccoon(Real hL, const MARScurvConfig<Order>& curvConfig,
                          vector<EdgeMark>& edgeMarks,
                          vector<vector<EdgeIndex>>& cyclesEdgesId,
                          vector<SmoothnessIndicator>& smoothConditions,
                          vector<vector<size_t>>& YinSetId);
};

OPTNONE_FUNC
inline auto
InterfaceGraphFactory::markRoseCurve(Point center, Real a, Real hL, rVec range,
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
OPTNONE_FUNC
auto InterfaceGraphFactory::markCurve(const Curve<DIM, Order>& crv, Real hL, Real lo, Real hi) {
  if (lo > hi) {
    lo = crv.getKnots().front();
    hi = crv.getKnots().back();
  }
  auto startpoint = crv(lo);
  auto endpoint = crv(hi);
  Real knotRange = hi - lo;
  auto eValue = [lo](auto& crv, Real t) { return crv(t + lo); };
  vector<Point> pts;
  vector<Real> ptsT;
  Real localHL = hL;
  Real step = localHL;
  const int sampleNum = 2;
  if (step > knotRange / sampleNum) step = knotRange / sampleNum;
  auto p0 = eValue(crv, 0);
  auto p1 = eValue(crv, step);
  double dist = norm(p0 - p1);
  while (dist < localHL / sampleNum) {
    step *= 2;
    p1 = eValue(crv, step);
    dist = norm(p0 - p1);
  }
  while (dist > 2 * localHL / sampleNum) {
    step /= 2;
    p1 = eValue(crv, step);
    dist = norm(p0 - p1);
  }
  int num = std::ceil(knotRange / step) * 2UL;
  step = knotRange / num;
  pts.push_back(startpoint);
  ptsT.push_back(0);
  for (int i = 0; i < num; i++) {
    auto pi1 = eValue(crv, step * (i + 1));
    Real remain = step;
    while (true) {
      if (norm(pi1 - pts.back()) < localHL) {
        pts.push_back(pi1);
        ptsT.push_back(step * (i + 1));
        break;
      }
      Real localStep = remain / 2;
      Real angle = step * (i + 1) - remain + localStep;
      Point local = eValue(crv, angle);
      while (norm(pts.back() - local) > localHL) {
        localStep /= 2;
        angle = step * (i + 1) - remain + localStep;
        local = eValue(crv, angle);
      }
      pts.push_back(local);
      remain -= localStep;
      ptsT.push_back(angle);
    }
  }
  if (norm(endpoint - pts.back()) < distTol()) {
    pts.pop_back();
    ptsT.pop_back();
  }
  pts.push_back(endpoint);
  ptsT.push_back(knotRange);

  vector<Point> ptsNew;
  vector<Real> ptsTNew;
  while (localHL > hL) {
    localHL /= 2;
    for (int i = 0; i < pts.size(); i++) {
      if (i == 0 || norm(pts[i] - pts[i - 1]) < localHL) {
        ptsNew.push_back(pts[i]);
        ptsTNew.push_back(ptsT[i]);
        continue;
      }
      Real localStep = (ptsT[i] - ptsT[i - 1]) / 2;
      Real angle = ptsT[i - 1];
      while (true) {
        auto localPt = eValue(crv, angle + localStep);
        while (norm(localPt - ptsNew.back()) > localHL) {
          localStep /= 2;
          localPt = eValue(crv, angle + localStep);
        }
        angle += localStep;
        ptsNew.push_back(localPt);
        ptsTNew.push_back(angle);
        localStep *= 2;
        if (norm(localPt - pts[i]) < localHL) {
          ptsNew.push_back(pts[i]);
          ptsTNew.push_back(ptsT[i]);
          break;
        }
      }
    }
    std::swap(pts, ptsNew);
    std::swap(ptsT, ptsTNew);
    ptsNew.clear();
    ptsTNew.clear();
  }

  return pts;
}

template <int Order>
OPTNONE_FUNC
auto InterfaceGraphFactory::markCurve(const Curve<DIM, Order>& crv, Real hL,
                                      Real t0, Real t1,
                                      const MARScurvConfig<Order-2>& curvConfig) {
  const auto& knots = crv.getKnots();
  if (t0 > t1) {
    t0 = knots.front();
    t1 = knots.back();
  }
  t0 = std::max(t0, knots.front());
  t1 = std::min(t1, knots.back());
  if (t0 >= t1) return std::vector<Real>{t0};

  std::vector<Real> params;
  params.reserve(2 * knots.size() + 2);
  params.push_back(t0);
  const int initNum = 5;
  double dt = (t1 - t0) / initNum;
  for (int k = 1; k < initNum; ++k) params.push_back(t0 + k * dt);
  // for (auto k : knots) {
  //   if (k > t0 && k < t1) {
  //     params.push_back(k);
  //   }
  // }
  params.push_back(t1);

  // std::vector<int> fixed;
  std::vector<int> fixed(params.size(), 0);
  fixed.front() = 1;
  fixed.back() = 1;

  auto distBounds = [&](const std::vector<Real>& curv,
                        std::vector<Real>& hLSeq, Real& lowerScale,
                        Real& upperScale, Real base) {
    if (curvConfig.used) {
      curvConfig.hlGenerator(curv, hLSeq, lowerScale, upperScale, base);
    } else {
      lowerScale = 0.2;
      upperScale = 0.6;
      hLSeq.assign(curv.size(), base);
    }
  };

  auto ts = crv.curvatureBoundParams(params, fixed, hL, distBounds);

  return ts;
  // std::vector<Point> pts;
  // pts.reserve(ts.size());
  // for (auto t : ts) pts.push_back(crv(t));
  // return pts;
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
OPTNONE_FUNC
auto InterfaceGraphFactory::createGraph41(Real hL)
    -> approxInterfaceGraph<Order> {
  vector<EdgeMark> edgeMarks(16);
  vector<SmoothnessIndicator> smoothConditions;
  vector<vector<EdgeIndex>> cyclesEdgesId(11);
  vector<vector<size_t>> YinSetId(6);
  // read tangle curve from file
  auto root = std::string(ROOT_DIR);
  std::ifstream infile(root + "/test/data1/Graph41.input");
  // auto readPts = [](auto& infile) {
  //   int num;
  //   infile.read((char*)&num, sizeof(num));
  //   vector<Real> xy(2UL * num);
  //   infile.read((char*)xy.data(), 2UL * num * sizeof(Real));
  //   vector<Point> pts;
  //   for (int i = 0; i < num; i++) {
  //     pts.push_back(Point{xy[2UL * i], xy[2UL * i + 1]});
  //   }
  //   return pts;
  // };

  auto crv1 = Curve<2, 6>::load(infile);
  auto crv2 = Curve<2, 6>::load(infile);
  Real k[4];
  infile.read((char*)k, sizeof(k) * 4UL);
  if (norm(crv1(k[0]) - crv2(k[3]), 2) > distTol() ||
      norm(crv1(k[1]) - crv2(k[2]), 2) > distTol()) {
    throw std::runtime_error("Invalid input file");
  }
  vector<Curve<2, 6>> splitRes;
  auto checkp0 = crv1(k[0]);
  k[0] = k[0] + crv1.getKnots().back() - crv1.getKnots().front();
  crv1.split({k[1]}, splitRes, distTol(), true);
  std::swap(crv1, splitRes[0]);
  if (Real v = norm(crv1(k[0]) - checkp0, 2); v > distTol()) {
    throw std::runtime_error("change k[0] fail.");
  }

  edgeMarks[0] = markCurve(crv1, hL, crv1.getKnots().front(), k[0]);
  edgeMarks[8] = markCurve(crv1, hL, k[0], crv1.getKnots().back());
  edgeMarks[12] = markCurve(crv2, hL, 0, k[2]);
  edgeMarks[11] = markCurve(crv2, hL, k[2], k[3]);
  edgeMarks[13] = markCurve(crv2, hL, k[3], crv2.getKnots().back());

  // Ellipse 1
  Point ell1Center{0.4, 0.575};
  rVec radius1{0.125, 0.09};
  edgeMarks[6] = markEllipse(ell1Center, radius1, hL, {0, M_PI / 2});
  edgeMarks[1] = markEllipse(ell1Center, radius1, hL, {M_PI / 2, M_PI});
  edgeMarks[3] = markEllipse(ell1Center, radius1, hL, {M_PI, 2 * M_PI});
  edgeMarks[5] = markSegment(
      Segment<2>(rVec{ell1Center[0] + radius1[0], ell1Center[1]}, ell1Center),
      hL);
  edgeMarks[2] = markSegment(
      Segment<2>(ell1Center, rVec{ell1Center[0] - radius1[0], ell1Center[1]}),
      hL);
  edgeMarks[4] = markSegment(
      Segment<2>(rVec{ell1Center[0], ell1Center[1] + radius1[1]}, ell1Center),
      hL);

  // Ellipse 2
  Point ell2Center{0.4, 0.415};
  rVec radius2{0.07, 0.035};
  edgeMarks[7] = markEllipse(ell2Center, radius2, hL, {0, 2 * M_PI});

  // Rose
  Point roseCenter{0.67, 0.475};
  Real radiusRose = 0.09;
  edgeMarks[14] = markRoseCurve(roseCenter, radiusRose, hL,
                                {M_PI * 4.0 / 6, M_PI * 5.0 / 6});
  std::reverse(edgeMarks[14].begin(), edgeMarks[14].end());
  edgeMarks[9] = markRoseCurve(roseCenter, radiusRose, hL,
                               {M_PI * 5.0 / 6, M_PI * 6.0 / 6});
  std::reverse(edgeMarks[9].begin(), edgeMarks[9].end());
  edgeMarks[10] = markRoseCurve(roseCenter, radiusRose, hL,
                                {M_PI * 6.0 / 6, M_PI * 7.0 / 6});
  std::reverse(edgeMarks[10].begin(), edgeMarks[10].end());
  edgeMarks[15] = markRoseCurve(roseCenter, radiusRose, hL,
                                {M_PI * 7.0 / 6, M_PI * 8.0 / 6});
  std::reverse(edgeMarks[15].begin(), edgeMarks[15].end());

  smoothConditions.emplace_back(13, 12);
  smoothConditions.emplace_back(12, 14);

  smoothConditions.emplace_back(1, 9);
  smoothConditions.emplace_back(9, 1);

  smoothConditions.emplace_back(7, 2);
  smoothConditions.emplace_back(2, 4);
  smoothConditions.emplace_back(4, 7);
  smoothConditions.emplace_back(6, 3);

  smoothConditions.emplace_back(8, 8);

  smoothConditions.emplace_back(16, 11);
  smoothConditions.emplace_back(11, 10);
  smoothConditions.emplace_back(10, 15);

  cyclesEdgesId[0] = (vector<EdgeIndex>{2, -3, -5});
  cyclesEdgesId[1] = (vector<EdgeIndex>{5, -6, 7});
  cyclesEdgesId[2] = (vector<EdgeIndex>{3, 4, 6});
  cyclesEdgesId[3] = (vector<EdgeIndex>{1, 9});
  cyclesEdgesId[4] = (vector<EdgeIndex>{-2, -7, -4});
  cyclesEdgesId[5] = (vector<EdgeIndex>{-8});
  cyclesEdgesId[6] = (vector<EdgeIndex>{13, 12, 14});
  cyclesEdgesId[7] = (vector<EdgeIndex>{-9, -12});
  cyclesEdgesId[8] = (vector<EdgeIndex>{10, 15});
  cyclesEdgesId[9] = (vector<EdgeIndex>{16, 11});
  cyclesEdgesId[10] = (vector<EdgeIndex>{-1, -13, -14});
  YinSetId[0] = (vector<size_t>{0});
  YinSetId[1] = (vector<size_t>{1});
  YinSetId[2] = (vector<size_t>{2});
  YinSetId[3] = (vector<size_t>{3, 4, 5, 6});
  YinSetId[4] = (vector<size_t>{7, 8, 9});
  YinSetId[5] = (vector<size_t>{10});

  return approxInterfaceGraph<Order>(std::move(edgeMarks), smoothConditions,
                                     std::move(cyclesEdgesId),
                                     std::move(YinSetId), distTol());
}

template <int Order>
OPTNONE_FUNC
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
  const Real T = params.field.T;
  const auto velocityPtr = params.field.velocity;

  // curvature-based adaption
  const bool curvUsed = params.curvature.curvUsed;
  const Real rCMin = params.curvature.rCMin;
  const Real rhoMin = params.curvature.rhoMin;
  const Real rhoMax = params.curvature.rhoMax;
  const auto sigma = [](Real x) { return x; };
  const InterfaceGraphFactory::MARScurvConfig<Order> curvConfig{
      .used = curvUsed,
      .rCMin = rCMin,
      .rhoCRange = Interval<1>{rhoMin, rhoMax},
      .sigma = sigma,
      .hlGenerator =
          HLGeneratorFactory::create(rCMin, rhoMin, rhoMax, rTiny, sigma),
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
  vector<Marsn2D::approxInterfaceGraph<Order>> vecDomain;
  vector<int> vecN;
  vector<Box<DIM>> vecBox;
  // vector<rVec> vecH;
  vector<Real> vecHL;
  vector<Real> vecDt;
  Real N = N0;
  auto hLComp = [&](Real N) {
    rVec h = (hi - lo) / N;
    Real hL = hLCoefficient * std::pow(min(h[0], h[1]), ((double)aimOrder) / Order);
    Real dt = std::min(h[0], h[1]) / uM * Cr;
    // Real initialDist = hL / 2;
    Real initialDist = hL / 2 * (curvUsed ? std::sqrt(rCMin) : 1);
    return make_tuple(h, hL, dt, initialDist);
  };
  for (uint i = 0; i < nGrid; ++i) {
    auto [h, hL, dt, initialDist] = hLComp(N);
    // rVec h = (hi - lo) / N;
    // Real hL = hLCoefficient * std::pow(min(h[0], h[1]), aimOrder / Order);
    // Real dt = std::min(h[0], h[1]) / uM * Cr;
    // Real initialDist = hL / 2 * (curvUsed ? std::sqrt(rCMin) : 1);
    if (params.domain.name == "Disk") {
      const auto disk = InterfaceGraphFactory::createDiskGraph<Order>(
          center, radius, initialDist, parts);
      vecDomain.emplace_back(std::move(disk));
    } else if (params.domain.name == "Graph41") {
      vecDomain.emplace_back(
          InterfaceGraphFactory::createGraph41<Order>(initialDist));
    } else if (params.domain.name == "Raccoon22") {
      string path = std::string(ROOT_DIR) + "/test/data1/Raccoon.input";
      // string path = std::string(ROOT_DIR) + "/test/data1/Raccoon.inputbak";
      vecDomain.emplace_back(InterfaceGraphFactory::inputFromSVG<Order>(
          path, hL / 2, curvConfig));
    } else if (params.domain.name == "Pig15") {
      string path = std::string(ROOT_DIR) + "/test/data1/Pig.input";
      vecDomain.emplace_back(InterfaceGraphFactory::inputFromSVG<Order>(
          path, hL / 2, curvConfig));
    }
    vecN.emplace_back(N);
    vecBox.emplace_back(0, N - 1);
    // vecH.emplace_back(h);
    vecHL.emplace_back(hL);
    vecDt.emplace_back(dt);
    N = N * 2;
  }

  // accurate solution
  auto [h, hL, dt, initialDist] = hLComp(N);
  // rVec h = (hi - lo) / (1 * N);
  // Real hL = hLCoefficient * std::pow(min(h[0], h[1]), aimOrder / Order);
  // Real initialDist = hL / 2 * (curvUsed ? std::sqrt(rCMin) : 1);
  auto exactDomain = InterfaceGraphFactory::createEmpty<Order>();
  if (params.domain.name == "Disk") {
    exactDomain = InterfaceGraphFactory::createDiskGraph<Order>(
        center, radius, initialDist, parts);
  } else if (params.domain.name == "Graph41") {
    exactDomain = InterfaceGraphFactory::createGraph41<Order>(initialDist);
  } else if (params.domain.name == "Raccoon22") {
    string path = std::string(ROOT_DIR) + "/test/data1/Raccoon.input";
    // string path = std::string(ROOT_DIR) + "/test/data1/Raccoon.inputbak";
    exactDomain = InterfaceGraphFactory::inputFromSVG<Order>(path, hL / 2,
                                                             curvConfig);
  } else if (params.domain.name == "Pig15") {
    string path = std::string(ROOT_DIR) + "/test/data1/Pig.input";
    exactDomain = InterfaceGraphFactory::inputFromSVG<Order>(path, hL / 2,
                                                             curvConfig);
  } else {
    throw runtime_error("unDeal shape.");
  }
  // MARSn2D<Order, VectorFunction> CM(timeIntegrator, hL, rTiny, curvConfig,
  //                                   false);
  // CM.trackInterface(
  //     *velocityPtr, exactDomain, t0, 0, t0,
  //     typename Marsn2D::MARSn2D<Order, VectorFunction>::PlotConfig());
  return make_tuple(vecDomain, exactDomain, radius, exactArea, exactLength,
                    vecBox, vecN, aimOrder, vecHL, rTiny, nGrid, curvConfig,
                    plotConfig, printDetail, t0, vecDt, te, T, timeIntegrator,
                    velocityPtr);
  }

  template <int Order>
  OPTNONE_FUNC inline auto InterfaceGraphFactory::inputFromSVG(
      const std::string& binPath, Real hL,
      const MARScurvConfig<Order>& curvConfig) -> approxInterfaceGraph<Order> {
    std::ifstream in(binPath, std::ios::binary);
    if (!in) {
      throw std::runtime_error("Failed to open file: " + binPath);
    }

    auto readBytes = [&](auto& val) {
      using T = std::decay_t<decltype(val)>;
      in.read(reinterpret_cast<char*>(&val), sizeof(T));
      if (!in) throw std::runtime_error("Unexpected EOF in " + binPath);
    };

    auto readVec = [&](auto& vec, size_t n) {
      using T = typename std::decay_t<decltype(vec)>::value_type;
      vec.resize(n);
      if (n == 0) return;
      in.read(reinterpret_cast<char*>(vec.data()), sizeof(T) * n);
      if (!in) throw std::runtime_error("Unexpected EOF in " + binPath);
    };

    // header
    char magic[4];
    in.read(magic, 4);
    if (std::string(magic, 4) != "GEOB") {
      throw std::runtime_error("Invalid magic in " + binPath);
    }
    uint32_t version;
    readBytes(version);
    if (version != 4) {
      throw std::runtime_error("Unsupported version in " + binPath);
    }

    uint32_t numV, numE, numCurves, numNB, numYin;
    readBytes(numV);
    readBytes(numE);
    readBytes(numCurves);
    readBytes(numNB);
    readBytes(numYin);

    // V
    std::vector<Real> vbuf(2 * numV);
    readVec(vbuf, vbuf.size());
    vector<Point> V(numV);
    for (uint32_t i = 0; i < numV; ++i) {
      V[i] = Point{vbuf[2 * i], vbuf[2 * i + 1]};
    }

    // Edges meta
    vector<SplineEdgeRef> edges(numE);
    for (uint32_t i = 0; i < numE; ++i) {
      int32_t v1, v2, idx, curveId;
      int8_t dir;
      Real t0, t1;
      readBytes(v1);
      readBytes(v2);
      readBytes(idx);
      readBytes(curveId);
      readBytes(dir);
      readBytes(t0);
      readBytes(t1);
      edges[i] = {v1, v2, idx, curveId, static_cast<int>(dir), t0, t1};
    }

    // Curves (Curve<2,Order>::load)
    vector<Curve<2, Order + 2>> curves(numCurves);
    for (uint32_t ci = 0; ci < numCurves; ++ci) {
      curves[ci] = Curve<2, Order + 2>::load(in);
    }

    // Smoothness
    vector<vector<SmoothPair>> smooth(numV);
    for (uint32_t vi = 0; vi < numV; ++vi) {
      uint32_t pairCount;
      readBytes(pairCount);
      if (pairCount == 0) continue;
      vector<int32_t> ebuf(2 * pairCount);
      vector<int8_t> fbuf(2 * pairCount);
      in.read(reinterpret_cast<char*>(ebuf.data()),
              sizeof(int32_t) * ebuf.size());
      in.read(reinterpret_cast<char*>(fbuf.data()),
              sizeof(int8_t) * fbuf.size());
      if (!in) throw std::runtime_error("Unexpected EOF in " + binPath);
      for (uint32_t k = 0; k < pairCount; ++k) {
        smooth[vi].push_back(SmoothPair{ebuf[2 * k + 0], fbuf[2 * k + 0],
                                        ebuf[2 * k + 1], fbuf[2 * k + 1]});
      }
    }

    // newBoundaries
    vector<BoundaryLoop> newBoundaries;
    for (uint32_t bi = 0; bi < numNB; ++bi) {
      uint8_t col[3];
      in.read(reinterpret_cast<char*>(col), 3);
      BoundaryLoop loop;
      loop.color = {static_cast<int>(col[0]), static_cast<int>(col[1]),
                    static_cast<int>(col[2])};
      uint32_t edgeCount;
      readBytes(edgeCount);
      loop.edges.resize(edgeCount);
      for (uint32_t ej = 0; ej < edgeCount; ++ej) {
        int32_t edgeId;
        int8_t dir;
        readBytes(edgeId);
        readBytes(dir);
        loop.edges[ej] = BoundaryEdge{edgeId, dir};
      }
      newBoundaries.push_back(std::move(loop));
    }
    // YinSet boundary索引
    vector<vector<size_t>> YinSetId;
    YinSetId.resize(numYin);
    for (uint32_t yi = 0; yi < numYin; ++yi) {
      uint8_t col[3];
      in.read(reinterpret_cast<char*>(col), 3);
      uint32_t idxCount;
      readBytes(idxCount);
      vector<uint32_t> tmp(idxCount);
      readVec(tmp, idxCount);
      YinSetId[yi] = vector<size_t>(tmp.begin(), tmp.end());
    }
    if (YinSetId.empty()) {
      YinSetId.push_back({});
    }

    // 采样样条生成 EdgeMark
    vector<EdgeMark> edgeMarks(numE);
    vector<vector<Real>> curveMarks(numCurves);
    for (int32_t ci = 0; ci < numCurves; ++ci) {
      curveMarks[ci] = markCurve(curves[ci], hL, 1, 0, curvConfig);
    }
    for (uint32_t ei = 0; ei < numE; ++ei) {
      const auto& ref = edges[ei];
      auto& pts = edgeMarks[ei];
      if (ref.curveId < 0 || ref.curveId >= (int)numCurves)
        throw std::runtime_error("inputFromSVG::invalid curveId.");
      if (ref.dir < 0 || ref.t0 > ref.t1)
        throw std::runtime_error("inputFromSVG::edge reverse");
      const auto& crv = curves[ref.curveId];
      const auto& crvMark = curveMarks[ref.curveId];
      auto iter0 = std::lower_bound(crvMark.begin(), crvMark.end(), ref.t0 + distTol());
      auto iter1 = std::lower_bound(iter0, crvMark.end(), ref.t1 - distTol());
      if (iter0 == crvMark.end() || iter1 == crvMark.end() || iter0 == iter1)
        throw std::runtime_error("inputFromSVG::fail to search t0, t1");
      pts.emplace_back(crv(ref.t0));
      while (iter0 != iter1) pts.emplace_back(crv(*iter0++));
      pts.emplace_back(crv(ref.t1));
    }

    // cyclesEdgesId 来自 newBoundaries (若无则 regions)
    vector<vector<EdgeIndex>> cyclesEdgesId(newBoundaries.size());
    for (size_t i = 0; i < newBoundaries.size(); ++i) {
      for (const auto& e : newBoundaries[i].edges) {
        cyclesEdgesId[i].push_back((int)(e.edgeId + 1) * (e.dir >= 0 ? 1 : -1));
      }
    }

    // smoothConditions (edge 下标 +1)
    vector<SmoothnessIndicator> smoothConditions;
    for (const auto& vlist : smooth) {
      for (const auto& p : vlist) {
        smoothConditions.emplace_back(p.eA + 1, p.eB + 1);
      }
    }

    auto name = binPath.substr(binPath.size()-13, 7);
    bool isRaccoon = name == "Raccoon";
    if (isRaccoon)
      dealRaccoon<Order>(hL, curvConfig, edgeMarks, cyclesEdgesId,
                         smoothConditions, YinSetId);

    return approxInterfaceGraph<Order>(
        std::move(edgeMarks), std::move(smoothConditions),
        std::move(cyclesEdgesId), std::move(YinSetId), distTol());
  }

template <int Order>
auto InterfaceGraphFactory::createEmpty()
    -> approxInterfaceGraph<Order> {
  vector<EdgeMark> edgeMarks;
  vector<SmoothnessIndicator> smoothConditions;
  vector<vector<EdgeIndex>> cyclesEdgesId;
  vector<vector<size_t>> YinSetId;
  return approxInterfaceGraph<Order>(
      std::move(edgeMarks), std::move(smoothConditions),
      std::move(cyclesEdgesId), std::move(YinSetId), distTol());
}

template <int Order>
auto InterfaceGraphFactory::dealRaccoon(
    Real hL, const MARScurvConfig<Order>& curvConfig,
    vector<EdgeMark>& edgeMarks, vector<vector<EdgeIndex>>& cyclesEdgesId,
    vector<SmoothnessIndicator>& smoothConditions,
    vector<vector<size_t>>& YinSetId) {
  Real initialDist = hL * (curvConfig.used ? std::sqrt(curvConfig.rCMin) : 1);
  Point roseCenter{0.45, 0.38};
  Real radiusRose = 0.06;
  int kEdge = edgeMarks.size();
  edgeMarks.resize(kEdge + 4);

  edgeMarks[kEdge] = markRoseCurve(roseCenter, radiusRose, initialDist,
                                   {M_PI * 4.0 / 6, M_PI * 5.0 / 6});
  std::reverse(edgeMarks[kEdge].begin(), edgeMarks[kEdge].end());
  edgeMarks[kEdge + 1] = markRoseCurve(roseCenter, radiusRose, initialDist,
                                       {M_PI * 5.0 / 6, M_PI * 6.0 / 6});
  std::reverse(edgeMarks[kEdge + 1].begin(), edgeMarks[kEdge + 1].end());
  edgeMarks[kEdge + 2] = markRoseCurve(roseCenter, radiusRose, initialDist,
                                       {M_PI * 6.0 / 6, M_PI * 7.0 / 6});
  std::reverse(edgeMarks[kEdge + 2].begin(), edgeMarks[kEdge + 2].end());
  edgeMarks[kEdge + 3] = markRoseCurve(roseCenter, radiusRose, initialDist,
                                       {M_PI * 7.0 / 6, M_PI * 8.0 / 6});
  std::reverse(edgeMarks[kEdge + 3].begin(), edgeMarks[kEdge + 3].end());

  smoothConditions.emplace_back(kEdge + 4, kEdge + 3);
  smoothConditions.emplace_back(kEdge + 3, kEdge + 2);
  smoothConditions.emplace_back(kEdge + 2, kEdge + 1);

  int kCyc = cyclesEdgesId.size();
  cyclesEdgesId.resize(kCyc + 2);
  cyclesEdgesId[kCyc] = (vector<EdgeIndex>{kEdge + 2, kEdge + 1});
  cyclesEdgesId[kCyc + 1] = (vector<EdgeIndex>{kEdge + 4, kEdge + 3});

  YinSetId[19].push_back(kCyc);
  YinSetId[19].push_back(kCyc + 1);
  // auto val = YinSetId[0];
  // YinSetId.clear();
  // YinSetId.push_back(val);
}

}  // namespace Marsn2D
