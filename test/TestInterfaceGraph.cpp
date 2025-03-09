#include "testHeader.h"

bool addInner = false;
bool output = false;

TEST_CASE("Disk 4 part interface.", "[InterfaceGraph][Disk][CutCell]") {
  Point center{0.5, 0.5};
  Real radio = 0.25;
  vector<Real> parts = {0.0, 0.25, 0.5, 0.75};
  for (auto& part : parts) part *= 2 * M_PI;
  // Real exactArea = M_PI * radio * radio / parts.size();
  // Real exactLength = 2 * M_PI * radio / parts.size() + 2 * radio;

  string name = "Circle";
  // Real lo = 0.0;
  // Real hi = 1.0;
  Real hL = 0.0001;

  auto dir = rootDir + "/results/InterfaceGraph/";
  mkdir(dir.c_str(), 0755);

  SECTION("2-th order disk") {
    constexpr int Order = 2;
    auto disk = Generator::createDisk<Order>(center, radio, hL, parts);
    int N = 4;
    Box<2> box(0, N - 1);
    // rVec h = (hi - lo) / N;
    auto yinSets = disk.approxYinSet();

    for (int i = 0; i < yinSets.size(); ++i) {
      string fileName = to_string(Order) + name + to_string(N) + "_" + to_string(i) + ".dat";
      std::ofstream of(dir + fileName, std::ios_base::binary);
      if (output) yinSets[i].dump(of);
    }
  }
}