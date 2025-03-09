#include "testHeader.h"

bool addInner = false;
bool output = false;

TEST_CASE("Circle Cut Cell", "[Circle][CutCell][SegmentedRealizableSpadjor]") {
  Point center{0.5, 0.5};
  Real radio = 0.25;
  Real exactArea = M_PI * radio * radio;
  Real exactLength = 2 * M_PI * radio;
  string name = "Circle";
  Real lo = 0.0;
  Real hi = 1.0;
  Real hL = 0.0001;
  Interval<2> range{lo, hi};

  SECTION("4-th order circle") {
    constexpr int Order = 4;
    auto circle = Generator::createEllipse<Order>(center, radio, hL);
    int N = 1024;
    Box<2> box(0, N - 1);
    rVec h = (hi - lo) / N;
    auto [res, boundary, tags] = circle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + name + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);

    THEN("Integral Area and length.") {
      Real totalArea = 0;
      Real length = 0;
      Real fullCell = h[0] * h[1];
      loop_box_2(res.box(), i0, i1) {
        if (tags(i0, i1) == 1) {
          totalArea += fullCell;
        } else if (tags(i0, i1) == 0) {
          for (const auto& crv : res(i0, i1)->getBoundaryCycles())
            totalArea += area(crv);
          for (const auto& crv : boundary(i0, i1)) length += arclength(crv);
        }
      }
      Real areaError = std::abs(totalArea - exactArea) / (exactArea);
      Real lengthError = std::abs(length - exactLength) / (exactLength);
      std::cout << name + " areaError = " << areaError << '\n';
      std::cout << name + " lengthError = " << lengthError << '\n';
      REQUIRE(std::fabs(areaError) < distTol() / exactArea);
      REQUIRE(std::fabs(lengthError) < distTol() / exactLength);
    }
  }

  SECTION("2-th order circle") {
    constexpr int Order = 2;
    auto circle = Generator::createEllipse<Order>(center, radio, 0.2);
    int N = 16;
    Box<2> box(0, N - 1);
    auto [res, boundary, tags] = circle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + "Circle" + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);
  }
}

TEST_CASE("Rectangle Cut Cell",
          "[Rectangle][CutCell][SegmentedRealizableSpadjor]") {
  Point rectLo = {0.25, 0.25};
  Point rectHi = {0.5, 0.75};
  Real exactArea = (rectHi[0] - rectLo[0]) * (rectHi[1] - rectLo[1]);
  Real exactLength = ((rectHi[0] - rectLo[0]) + (rectHi[1] - rectLo[1])) * 2;
  string name = "Rectangle";
  Real theta = 0;
  bool orientation = true;
  Real hL = 0.001;
  Real lo = 0.0;
  Real hi = 1.0;

  SECTION("4-th order rectangle") {
    constexpr int Order = 4;
    CurveFactory<2, Order> crvFactory;
    Interval<2> range{lo, hi};
    int N = 512;
    rVec h = (hi - lo) / N;
    Box<2> box(0, N - 1);
    stringstream is;
    is << name << ' ';
    is << rectLo[0] << ' ' << rectLo[1] << ' ' << rectHi[0] << ' ' << rectHi[1]
       << ' ' << theta << ' ' << orientation << ' ' << hL;
    auto rect = *crvFactory.createCurve(is.str());
    vector<OrientedJordanCurve<2, Order>> curves({rect});
    YinSet<DIM, Order> rectangle({curves}, distTol());
    auto [res, boundary, tags] = rectangle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + name + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);

    THEN("Integral Area and length.") {
      Real totalArea = 0;
      Real length = 0;
      Real fullCell = h[0] * h[1];
      loop_box_2(res.box(), i0, i1) {
        if (tags(i0, i1) == 1) {
          totalArea += fullCell;
        } else if (tags(i0, i1) == 0) {
          for (const auto& crv : res(i0, i1)->getBoundaryCycles())
            totalArea += area(crv);
          for (const auto& crv : boundary(i0, i1)) length += arclength(crv);
        }
      }
      Real areaError = std::abs(totalArea - exactArea) / (exactArea);
      Real lengthError = std::abs(length - exactLength) / (exactLength);
      std::cout << name + " areaError = " << areaError << '\n';
      std::cout << name + " lengthError = " << lengthError << '\n';
      REQUIRE(std::fabs(areaError) < distTol() / exactArea);
      REQUIRE(std::fabs(lengthError) < distTol() / exactLength);
    }
  }

  SECTION("2-th order rectangle") {
    constexpr int Order = 2;
    CurveFactory<2, Order> crvFactory;
    Interval<2> range{lo, hi};
    int N = 512;
    Box<2> box(0, N - 1);
    stringstream is;
    is << name << ' ';
    is << rectLo[0] << ' ' << rectLo[1] << ' ' << rectHi[0] << ' ' << rectHi[1]
       << ' ' << theta << ' ' << orientation << ' ' << hL;
    auto rect = *crvFactory.createCurve(is.str());
    vector<OrientedJordanCurve<2, Order>> curves({rect});
    YinSet<DIM, Order> rectangle({curves}, distTol());
    auto [res, boundary, tags] = rectangle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + name + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);
  }
}