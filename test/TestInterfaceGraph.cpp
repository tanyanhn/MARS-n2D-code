#include "testHeader.h"

bool addInner = false;
bool output = true;

TEST_CASE("Disk 4 part interface.", "[InterfaceGraph][Disk][CutCell][4]") {
  Point center{0.5, 0.5};
  Real radio = 0.25;
  vector<Real> parts = {0.0, 0.25, 0.5, 0.75};
  for (auto& part : parts) part *= 2 * M_PI;
  Real exactArea = M_PI * radio * radio / parts.size();
  Real exactLength = 2 * M_PI * radio / parts.size() + 2 * radio;

  string name = "Circle";
  Real lo = 0.0;
  Real hi = 1.0;
  Real hL = 0.0001;
  Interval<2> range{lo, hi};

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
      auto [res, boundary, tags] = yinSets[i].cutCell(box, range, addInner);
      if (output) {
        string fileName = to_string(Order) + name + to_string(N) + "_" +
                          to_string(i + 1) + ".dat";
        std::ofstream of(dir + fileName, std::ios_base::binary);
        dumpTensorYinSet<Order>(res, of);
        // yinSets[i].dump(of);
      }
    }
  }

  SECTION("4-th order disk") {
    constexpr int Order = 4;
    auto disk = Generator::createDisk<Order>(center, radio, hL, parts);
    int N = 4;
    Box<2> box(0, N - 1);
    rVec h = (hi - lo) / N;
    auto yinSets = disk.approxYinSet();

    for (int i = 0; i < yinSets.size(); ++i) {
      auto [res, boundary, tags] = yinSets[i].cutCell(box, range, addInner);
      if (output) {
        string fileName = to_string(Order) + name + to_string(N) + "_" +
                          to_string(i + 1) + ".dat";
        std::ofstream of(dir + fileName, std::ios_base::binary);
        dumpTensorYinSet<Order>(res, of);
        // yinSets[i].dump(of);
      }

      INFO("Integral Area and length.");
      Real totalArea = 0;
      Real length = 0;
      Real fullCell = h[0] * h[1];
      loop_box_2(res.box(), i0, i1) {
        if (tags(i0, i1) == 1) {
          totalArea += fullCell;
        } else if (tags(i0, i1) == 0) {
          if (res(i0, i1)) {
            for (const auto& crv : res(i0, i1)->getBoundaryCycles())
              totalArea += area(crv);
          }
          for (const auto& crv : boundary(i0, i1)) length += arclength(crv);
        }
      }
      if (i == yinSets.size() - 1) {
        exactArea = 1 - exactArea * (yinSets.size() - 1);
        exactLength = 2 * M_PI * radio;
      }
      Real areaError = std::abs(totalArea - exactArea) / (exactArea);
      Real lengthError = std::abs(length - exactLength) / (exactLength);
      std::cout << name + to_string(i) + " areaError = " << areaError << '\n';
      std::cout << name + to_string(i) + " lengthError = " << lengthError
                << '\n';
      REQUIRE(std::fabs(areaError) < distTol() / exactArea);
      REQUIRE(std::fabs(lengthError) < distTol() / exactLength);
    }
  }
}

TEST_CASE("Disk 5 part interface.", "[InterfaceGraph][Disk][CutCell][5]") {
  Point center{0.5, 0.5};
  Real radio = 0.25;
  vector<Real> parts = {1.0 / 20, 5.0 / 20, 9.0 / 20, 13.0 / 20, 17.0 / 20};
  for (auto& part : parts) part *= 2 * M_PI;
  Real exactArea = M_PI * radio * radio / parts.size();
  Real exactLength = 2 * M_PI * radio / parts.size() + 2 * radio;

  string name = "Circle";
  Real lo = 0.0;
  Real hi = 1.0;
  Real hL = 0.0001;
  Interval<2> range{lo, hi};

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
      auto [res, boundary, tags] = yinSets[i].cutCell(box, range, addInner);
      if (output) {
        string fileName = to_string(Order) + name + to_string(N) + "_" +
                          to_string(i + 1) + ".dat";
        std::ofstream of(dir + fileName, std::ios_base::binary);
        dumpTensorYinSet<Order>(res, of);
        // yinSets[i].dump(of);
      }
    }
  }

  SECTION("4-th order disk") {
    constexpr int Order = 4;
    auto disk = Generator::createDisk<Order>(center, radio, hL, parts);
    int N = 4;
    Box<2> box(0, N - 1);
    rVec h = (hi - lo) / N;
    auto yinSets = disk.approxYinSet();

    for (int i = 1; i < yinSets.size(); ++i) {
      auto [res, boundary, tags] = yinSets[i].cutCell(box, range, addInner);
      if (output) {
        string fileName = to_string(Order) + name + to_string(N) + "_" +
                          to_string(i + 1) + ".dat";
        std::ofstream of(dir + fileName, std::ios_base::binary);
        dumpTensorYinSet<Order>(res, of);
        // yinSets[i].dump(of);
      }

      INFO("Integral Area and length.");
      Real totalArea = 0;
      Real length = 0;
      Real fullCell = h[0] * h[1];
      loop_box_2(res.box(), i0, i1) {
        if (tags(i0, i1) == 1) {
          totalArea += fullCell;
        } else if (tags(i0, i1) == 0) {
          if (res(i0, i1)) {
            if (res(i0, i1)) {
              for (const auto& crv : res(i0, i1)->getBoundaryCycles())
                totalArea += area(crv);
            }
          }
          for (const auto& crv : boundary(i0, i1)) length += arclength(crv);
        }
      }
      if (i == yinSets.size() - 1) {
        exactArea = 1 - exactArea * (yinSets.size() - 1);
        exactLength = 2 * M_PI * radio;
      }
      Real areaError = std::abs(totalArea - exactArea) / (exactArea);
      Real lengthError = std::abs(length - exactLength) / (exactLength);
      std::cout << name + to_string(i) + " areaError = " << areaError <<
      '\n'; std::cout << name + to_string(i) + " lengthError = " <<
      lengthError
                << '\n';
      REQUIRE(std::fabs(areaError) < distTol() / exactArea);
      REQUIRE(std::fabs(lengthError) < distTol() / exactLength);
    }
  }
}