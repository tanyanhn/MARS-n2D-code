// Only include by test file.
#include "Core/Config.h"
#include "Core/Vec.h"
#include "Core/dirConfig.h"
#include "Marsn2D/InterfaceGraph.h"
#include "YinSet/YinSet.h"
#include "catch_amalgamated.hpp"

using namespace std;
using namespace Catch;
using namespace MARSn2D;

using Point = Vec<Real, 2>;
using rVec = Point;
template <int Order>
using YinSetPtr = YinSet<2, Order>::YinSetPtr;

const std::string rootDir(ROOT_DIR);

struct Generator {
  static VecCompare<Real, 2> vCmp;
  template <int l, int r>
  [[nodiscard]] static auto randomCreateReal() -> Real;

  template <int Order>
  static auto createEllipse(Point center, rVec radio, Real hL,
                            rVec range = {0, 2 * M_PI}) -> YinSet<2, Order>;

  template <int Order>
  static auto createDiskGraph(Point center, rVec radio, Real hL, vector<Real> parts)
      -> approxInterfaceGraph<Order>;
  static auto markEllipse(Point center, rVec radio, Real hL, rVec range);
  static auto markSegment(Segment<DIM> seg, Real hL);
};

const Generator generator;
VecCompare<Real, 2> Generator::vCmp(distTol());

template <int l, int r>
auto Generator::randomCreateReal() -> Real {
  STATIC_CHECK(l < r);
  return GENERATE(take(1, random(l, r)));
}


inline auto Generator::markEllipse(Point center, rVec radius, Real hL, rVec range) {
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
                           vector<Real> parts) -> approxInterfaceGraph<Order> {
  using SmoothnessIndicator = InterfaceGraph::SmoothnessIndicator;
  using EdgeIndex = InterfaceGraph::EdgeIndex;
  using std::vector;

  vector<EdgeMark> edgeMarks;
  vector<SmoothnessIndicator> smoothConditions;
  vector<vector<EdgeIndex>> cyclesEdgesId(parts.size() + 1);
  vector<vector<size_t>> YinSetId(parts.size() + 1);

  auto eValue = [center, radius](Real angle) {
    return Point{center + rVec{radius[0] * cos(angle), radius[1] * sin(angle)}};
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

  return approxInterfaceGraph<Order>(std::move(edgeMarks), smoothConditions,
                                     std::move(cyclesEdgesId),
                                     std::move(YinSetId), distTol());
}

template <int Order>
void dumpTensorYinSet(const Tensor<YinSetPtr<Order>, 2>& data,
                      std::ostream& os) {
  int num = 0;
  loop_box_2(data.box(), i0, i1) {
    if (data(i0, i1)) {
      num++;
    }
  }
  os.write((char*)&num, sizeof(num));
  loop_box_2(data.box(), i0, i1) {
    if (data(i0, i1)) {
      data(i0, i1)->dump(os);
    }
  }
}