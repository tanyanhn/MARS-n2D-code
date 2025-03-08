#include "testHeader.h"

bool addInner = false;
bool output = false;

TEST_CASE("Circle Cut Cell", "[Circle][CutCell][SegmentedRealizableSpadjor]") {
  Point center{0.5, 0.5};
  Real radio = 0.3;

  SECTION("4-th order circle") {
    constexpr int Order = 4;
    auto circle = Generator::createEllipse<Order>(center, radio, 0.00001);
    int N = 1024;
    Real lo = 0.0;
    Real hi = 1.0;
    Interval<2> range{lo, hi};
    Box<2> box(0, N - 1);
    auto [res, tags] = circle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + "Circle" + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);
  }
  
  SECTION("2-th order circle") {
    constexpr int Order = 2;
    auto circle = Generator::createEllipse<Order>(center, radio, 0.2);
    int N = 16;
    Real lo = 0.0;
    Real hi = 1.0;
    Interval<2> range{lo, hi};
    Box<2> box(0, N - 1);
    auto [res, tags] = circle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + "Circle" + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);
  }
}


TEST_CASE("Rectangle Cut Cell", "[Rectangle][CutCell][SegmentedRealizableSpadjor]") {
  Point rectLo = {0.25, 0.25};
  Point rectHi = {0.5, 0.75};
  string name = "Rectangle";
  Real theta = 0;
  bool orientation = true;
  Real hL = 0.00001;
  Real lo = 0.0;
  Real hi = 1.0;

  SECTION("4-th order rectangle") {
    constexpr int Order = 4;
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
    auto [res, tags] = rectangle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + name + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);
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
    auto [res, tags] = rectangle.cutCell(box, range, addInner);
    // output res
    auto dir = rootDir + "/results/CutCell/";
    mkdir(dir.c_str(), 0755);
    string fileName = to_string(Order) + name + to_string(N) + ".dat";
    std::ofstream of(dir + fileName, std::ios_base::binary);
    if (output) dumpTensorYinSet<Order>(res, of);
  }
}