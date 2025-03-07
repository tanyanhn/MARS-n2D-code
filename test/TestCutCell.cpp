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