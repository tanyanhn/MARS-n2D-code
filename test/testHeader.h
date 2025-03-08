// Only include by test file.
#include "Core/Config.h"
#include "Core/Vec.h"
#include "Core/dirConfig.h"
#include "YinSet/YinSet.h"
#include "catch_amalgamated.hpp"

using namespace std;
using namespace Catch;

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
  static auto createEllipse(Point center, rVec radio, Real hL)
      -> YinSet<2, Order>;
};

const Generator generator;
VecCompare<Real, 2> Generator::vCmp(distTol());

template <int l, int r>
auto Generator::randomCreateReal() -> Real {
  STATIC_CHECK(l < r);
  return GENERATE(take(1, random(l, r)));
}

template <int Order>
auto Generator::createEllipse(Point center, rVec radio, Real hL)
    -> YinSet<2, Order> {
  vector<Point> pts;
  Real stepAngle = hL / (radio[0] + radio[1]) * 2;
  int step = 2 * M_PI / stepAngle;
  pts.reserve(step + 1);
  for (int i = 0; i < step; i++) {
    pts.push_back({center[0] + radio[0] * cos(stepAngle * i),
                   center[1] + radio[1] * sin(stepAngle * i)});
  }
  if (vCmp(pts.front(), pts.back()) != 0) pts.push_back(pts.front());

  auto crv = fitCurve<Order>(pts, Curve<2, Order>::periodic);
  vector<OrientedJordanCurve<2, Order>> vCrv{crv};
  YinSet<2, Order> YS(SegmentedRealizableSpadjor<Order>(vCrv), distTol());
  return YS;
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