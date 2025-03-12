#include <cmath>

#include "testHeader.h"

bool addInner = true;
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
    auto disk = Generator::createDiskGraph<Order>(center, radio, hL, parts);
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
    auto disk = Generator::createDiskGraph<Order>(center, radio, hL, parts);
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
      std::cout << name + to_string(i) + " areaError = " << areaError <<
      '\n'; std::cout << name + to_string(i) + " lengthError = " <<
      lengthError
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
    auto disk = Generator::createDiskGraph<Order>(center, radio, hL, parts);
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
    auto disk = Generator::createDiskGraph<Order>(center, radio, hL, parts);
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

TEST_CASE("Ellipse single.", "[InterfaceGraph][Ellipse][CutCell][1]") {
  Point center{0.5, 0.5};
  rVec radius = {0.3, 0.2};

  string name = "Ellipse";
  Real lo = 0.0;
  Real hi = 1.0;
  Real hL = 0.00005;
  Interval<2> range{lo, hi};

  auto dir = rootDir + "/results/InterfaceGraph/";
  mkdir(dir.c_str(), 0755);

  SECTION("2-th order ellipse") {
    constexpr int Order = 2;
    auto yinSet = Generator::createEllipse<Order>(center, radius, hL);
    std::vector<YinSet<DIM, Order>> yinSets;
    yinSets.push_back(yinSet);
    int N = 16;
    Box<2> box(0, N - 1);

    //=============================================
    // check the area
    INFO("Check integral area.");
    auto [res, boundary, tags] = yinSets[0].cutCell(box, range, addInner);
    Real exactArea = M_PI * radius[0] * radius[1];
    Real totalArea = 0;
    rVec h = (hi - lo) / N;
    Real fullCell = h[0] * h[1];
    loop_box_2(res.box(), i0, i1) {
      if (tags(i0, i1) == 1) {
        totalArea += fullCell;
      } else if (tags(i0, i1) == 0) {
        if (res(i0, i1)) {
          for (const auto& crv : res(i0, i1)->getBoundaryCycles())
            totalArea += area(crv);
        }
      }
    }
    Real areaError = std::abs(totalArea - exactArea) / (exactArea);
    std::cout << name + " areaError = " << areaError << '\n';
    //=============================================
    // dump
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
  }  // end of 2-th order

  SECTION("4-th order ellipse") {
    constexpr int Order = 4;
    auto yinSet = Generator::createEllipse<Order>(center, radius, hL);
    std::vector<YinSet<DIM, Order>> yinSets;
    yinSets.push_back(yinSet);
    int N = 16;
    Box<2> box(0, N - 1);

    //=============================================
    // check the area
    INFO("Check integral area.");
    auto [res, boundary, tags] = yinSets[0].cutCell(box, range, addInner);
    Real exactArea = M_PI * radius[0] * radius[1];
    Real totalArea = 0;
    rVec h = (hi - lo) / N;
    Real fullCell = h[0] * h[1];
    loop_box_2(res.box(), i0, i1) {
      if (tags(i0, i1) == 1) {
        totalArea += fullCell;
      } else if (tags(i0, i1) == 0) {
        if (res(i0, i1)) {
          for (const auto& crv : res(i0, i1)->getBoundaryCycles())
            totalArea += area(crv);
        }
      }
    }
    Real areaError = std::abs(totalArea - exactArea) / (exactArea);
    std::cout << name + " areaError = " << areaError << '\n';
    //=============================================
    // dump
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
  }  // end of 4-th order
}

TEST_CASE("Ellipse 3 part interface.",
          "[InterfaceGraph][Ellipse][CutCell][3]") {
  Point center{0.5, 0.5};
  rVec radius = {0.3, 0.2};
  vector<Real> parts = {0.0, 0.25, 0.5};
  for (auto& part : parts) part *= 2 * M_PI;
  // calculate the area
  Real ellipseArea = M_PI * radius[0] * radius[1];
  std::function<Real(Real, Real)> calculate_area = [&](Real theta1, Real
  theta2) -> Real {
    assert(theta2 > theta1);
    if(abs(theta2 - 2*M_PI) < 1e-10){
      return ellipseArea - calculate_area(0, theta1);
    }
    auto a = radius[0];
    auto b = radius[1];
    auto theta1_conv = atan2(a * sin(theta1), b * cos(theta1));
    auto theta2_conv = atan2(a * sin(theta2), b * cos(theta2));
    auto result = a * b / 2.0 * (theta2_conv - theta1_conv);
    return result;
  };
  vector<Real> exactAreas(parts.size());
  for (size_t i = 0; i < parts.size(); ++i) {
    if (i == parts.size() - 1) {
      exactAreas.at(i) = calculate_area(parts.at(i), 2 * M_PI);
    } else {
      exactAreas.at(i) = calculate_area(parts.at(i), parts.at(i + 1));
    }
  }

  string name = "Ellipse";
  Real lo = 0.0;
  Real hi = 1.0;
  Real hL = 0.0001;
  Interval<2> range{lo, hi};

  auto dir = rootDir + "/results/InterfaceGraph/";
  mkdir(dir.c_str(), 0755);

  SECTION("2-th order Ellipse") {
    constexpr int Order = 2;
    auto ellipse = Generator::createDiskGraph<Order>(center, radius, hL,
    parts); int N = 4; Box<2> box(0, N - 1);
    // rVec h = (hi - lo) / N;
    auto yinSets = ellipse.approxYinSet();

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
  } // end of 2-th order ellipse

  SECTION("4-th order ellipse") {
    constexpr int Order = 4;
    auto ellipse = Generator::createDiskGraph<Order>(center, radius, hL,
    parts); int N = 8; Box<2> box(0, N - 1); rVec h = (hi - lo) / N; auto
    yinSets = ellipse.approxYinSet();

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
      Real fullCell = h[0] * h[1];
      loop_box_2(res.box(), i0, i1) {
        if (tags(i0, i1) == 1) {
          totalArea += fullCell;
        } else if (tags(i0, i1) == 0) {
          if (res(i0, i1)) {
            for (const auto& crv : res(i0, i1)->getBoundaryCycles())
              totalArea += area(crv);
          }
        }
      }
      Real exactArea;
      if (i == yinSets.size() - 1) {
        exactArea = 1 - ellipseArea;
      }else{
        exactArea = exactAreas.at(i);
      }
      Real areaError = std::abs(totalArea - exactArea) / (exactArea);
      // Real lengthError = std::abs(length - exactLength) / (exactLength);
      std::cout << name + to_string(i) + " areaError = " << areaError <<
      '\n';
      // std::cout << name + to_string(i) + " lengthError = " << lengthError
                // << '\n';
      REQUIRE(std::fabs(areaError) < distTol() / exactArea);
      // REQUIRE(std::fabs(lengthError) < distTol() / exactLength);
    }
  } // end of 4-th order loop
}

TEST_CASE("Rose Curve.", "[InterfaceGraph][RoseCurve][CutCell][1]") {
  Point center{0.5, 0.5};
  Real a = 0.25;
  vector<Real> parts = {0.0, 1.0 / 3, 2.0 / 3};
  for (auto& part : parts) part *= M_PI;
  Real exactArea = M_PI * a * a / 12;
  // Real exactLength = 2 * M_PI * radio / parts.size() + 2 * radio;

  string name = "Rose";
  Real lo = 0.0;
  Real hi = 1.0;
  Real hL = 0.0005;
  Interval<2> range{lo, hi};

  auto dir = rootDir + "/results/InterfaceGraph/";
  mkdir(dir.c_str(), 0755);

  SECTION("2-th order rose") {
    constexpr int Order = 2;
    auto disk = Generator::createRoseGraph<Order>(center, a, hL, parts);
    int N = 16;
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

  SECTION("4-th order rose") {
    constexpr int Order = 4;
    auto disk = Generator::createRoseGraph<Order>(center, a, hL, parts);
    int N = 16;
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
      // Real length = 0;
      Real fullCell = h[0] * h[1];
      loop_box_2(res.box(), i0, i1) {
        if (tags(i0, i1) == 1) {
          totalArea += fullCell;
        } else if (tags(i0, i1) == 0) {
          if (res(i0, i1)) {
            for (const auto& crv : res(i0, i1)->getBoundaryCycles())
              totalArea += area(crv);
          }
          // for (const auto& crv : boundary(i0, i1)) length += arclength(crv);
        }
      }
      if (i == yinSets.size() - 1) {
        exactArea = 1 - exactArea * (yinSets.size() - 1);
        // exactLength = 2 * M_PI * radio;
      }
      Real areaError = std::abs(totalArea - exactArea) / (exactArea);
      // Real lengthError = std::abs(length - exactLength) / (exactLength);
      std::cout << name + to_string(i) + " areaError = " << areaError << '\n';
      // std::cout << name + to_string(i) + " lengthError = " <<
      // lengthError
      //           << '\n';
      REQUIRE(std::fabs(areaError) < distTol() / exactArea);
      // REQUIRE(std::fabs(lengthError) < distTol() / exactLength);
    }
  }
}