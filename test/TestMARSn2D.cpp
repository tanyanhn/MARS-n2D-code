#include <nlohmann/json.hpp>

#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/ERK.h"
#include "InterfaceTracking/RKButcher.h"
#include "InterfaceTracking/VelocityField.h"
#include "Marsn2D/InterfaceGraph.h"
#include "Marsn2D/MARSReadJson.hpp"
#include "Marsn2D/MARSn2D.h"
#include "Recorder/Timer.h"
#include "testHeader.h"
using namespace nlohmann;

bool addInner = false;
bool output = false;

// auto read_config(const std::string& filename) {
//   std::ifstream ifs(filename);
//   return json::parse(ifs);
// }

template <int Order>
auto diskTEST(const std::string& jsonFile) {
  auto params = read_config(jsonFile);
  distTol(params.tol);
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
  // const rVec h0 = (hi - lo) / N0;
  // const Real hL0 =
  //     hLCoefficient * std::pow(h0[0] * h0[1], 0.5 * aimOrder / Order);
  // const Box<DIM> box(0, N0 - 1);
  // const Real dt0 = std::min(h0[0], h0[1]) / uM * Cr;
  std::shared_ptr<TimeIntegrator<DIM, VectorFunction>> timeIntegrator;
  if (Order == 4) {
    timeIntegrator = make_shared<ERK<DIM, RK::ClassicRK4, VectorFunction>>();
  } else if (Order == 6) {
    timeIntegrator = make_shared<ERK<DIM, RK::Verner6, VectorFunction>>();
  } else if (Order == 8) {
    timeIntegrator =
        make_shared<ERK<DIM, RK::PrinceDormand8, VectorFunction>>();
  }
  const Real T = te;
  // const auto vortex = std::make_shared<Vortex>(T);
  // const auto velocityPtr = vortex;
  const auto velocityPtr = params.field.velocity;

  // curvature-based adaption
  const bool curvUsed = params.curvature.curvUsed;
  const Real rCMin = params.curvature.rCMin;
  const Real rhoMin = params.curvature.rhoMin;
  const Real rhoMax = params.curvature.rhoMax;
  const typename MARSn2D<Order, VectorFunction>::CurvatureAdaptionConfig
      curvConfig{
          .used = curvUsed,
          .rCMin = rCMin,
          .rhoCRange = Interval<1>{rhoMin, rhoMax},
          .sigma = [](Real x) { return x; },
      };

  // plot
  const int plotOutput = params.plotConfig.plotOutput;
  const bool plotInner = params.plotConfig.plotInner;
  const Real plotStride = params.plotConfig.plotStride;
  const bool printDetail = params.plotConfig.printDetail;
  const struct MARSn2D<Order, VectorFunction>::PlotConfig plotConfig {
    .output = plotOutput, .plotInner = plotInner,
    .opStride = (int)std::round(plotStride * te * uM / Cr),
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
  vector<approxInterfaceGraph<Order>> vecDisk;
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
    const auto disk =
        Generator::createDiskGraph<Order>(center, radius, hL, parts);
    vecDisk.emplace_back(std::move(disk));
    vecN.emplace_back(N);
    vecBox.emplace_back(0, N - 1);
    // vecH.emplace_back(h);
    vecHL.emplace_back(hL);
    vecDt.emplace_back(dt);
    N = N * 2;
  }

  // accurate solution
  N *= 2;
  rVec h = (hi - lo) / N;
  Real hL = hLCoefficient * std::pow(h[0] * h[1], 0.5 * aimOrder / Order);
  auto exactDisk = Generator::createDiskGraph<Order>(center, radius, hL, parts);

  return make_tuple(vecDisk, exactDisk, radius, exactArea, exactLength, vecBox,
                    vecN, aimOrder, vecHL, rTiny, nGrid, curvConfig, plotConfig,
                    printDetail, t0, vecDt, te, timeIntegrator, velocityPtr);
  // }
}

template <int Order>
void checkResult(auto& yinSets, auto& box, auto& range, auto& addInner,
                 auto& radius, auto& N, auto& h, auto& output, auto& dir,
                 auto& name, auto& exactLength, auto& exactArea) {
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
      if (tags[i0][i1] == 1 || tags[i0][i1] == 2) {
        totalArea += fullCell;
      } else if (tags[i0][i1] == 0 || tags[i0][i1] == -2 || tags[i0][i1] == 2) {
        if (res(i0, i1)) {
          for (const auto& crv : res(i0, i1)->getBoundaryCycles())
            totalArea += area(crv);
        }
        for (const auto& crv : boundary(i0, i1)) length += arclength(crv);
      }
    }
    if (i == yinSets.size() - 1) {
      exactArea = 1 - exactArea * (yinSets.size() - 1);
      exactLength = 2 * M_PI * radius[0];
    }
    Real areaError = std::abs(totalArea - exactArea) / (exactArea);
    Real lengthError = std::abs(length - exactLength) / (exactLength);
    std::cout << to_string(i) + "th areaError = " << areaError << '\n';
    std::cout << to_string(i) + "th lengthError = " << lengthError << '\n';
    // REQUIRE(std::fabs(areaError) < distTol() / exactArea);
    // REQUIRE(std::fabs(lengthError) < distTol() / exactLength);
  }
}

// TEST_CASE("Disk 4 vortex4, RK4: Area and Length Test.",
//           "[Disk][Vortex][4][MARSn2D]") {
//   ::Timer t("disk4Vortex4");
//   const auto* testName = "Disk4Vortex4";
//   auto dir = rootDir + "/results/TrackInterface/" + testName + "/";
//   mkdir(dir.c_str(), 0755);
//   SECTION("2-th order disk") {
//     constexpr int Order = 2;
//     auto [disk, radius, exactArea, exactLength, box, N, aimOrder, h, hL,
//     rTiny,
//           nGrid, curvConfig, plotConfig, printDetail, t0, dt, te,
//           timeIntegrator, vortex] =
//         diskTEST<Order>(rootDir + "/test/config/" + testName + ".json");
//     plotConfig.fName = dir + plotConfig.fName;
//     MARSn2D<Order, VectorFunction> CM(timeIntegrator, hL, rTiny, curvConfig,
//                                       printDetail);

//     CM.trackInterface(vortex, disk, t0, dt, te, plotConfig);

//     INFO("Integral Area and length.");
//     auto yinSets = disk.approxYinSet();

//     checkResult<Order>(yinSets, box, plotConfig.range, addInner, radius, N,
//     h,
//                        output, dir, plotConfig.fName, exactLength,
//                        exactArea);
//     t.~Timer();
//     ::Timer::printStatistics();
//   }
// }

TEST_CASE("Disk 4 vortex T = 4/8/12/16, Convergence Test.",
          "[Disk][Vortex][MARSn2D][Convergence]") {
  // const auto* testName = "Disk4Vortex4";
  const auto* testName = "Disk4Vortex8";
  // const auto* testName = "Disk4Vortex12";
  // const auto* testName = "Disk4Vortex16";
  ::Timer t(testName);
  auto dir = rootDir + "/results/TrackInterface/" + testName + "/";
  mkdir(dir.c_str(), 0755);
  RecorderInitialize(RecorderInfo{DebugLevel::INFO, dir});

  SECTION("4-th order disk") {
    constexpr int Order = 4;
    pushLogStage("Initialize");
    auto [vecDisk, exactDisk, radius, exactArea, exactLength, vecBox, vecN,
          aimOrder, vecHL, rTiny, nGrid, curvConfig, plotConfig, printDetail,
          t0, vecDt, te, timeIntegrator, velocityPtr] =
        diskTEST<Order>(rootDir + "/test/config/" + testName + ".json");
    popLogStage();

    pushLogStage("TrackInterface");
    for (uint i = 0; i < nGrid; ++i) {
      plotConfig.fName = to_string(Order) + "Circle" + "_grid" +
                         to_string(vecN[i]) + "_T" + to_string((int)te),
      plotConfig.fName = dir + plotConfig.fName;
      MARSn2D<Order, VectorFunction> CM(timeIntegrator, vecHL[i], rTiny,
                                        curvConfig, printDetail);

      CM.trackInterface(*velocityPtr, vecDisk[i], t0, vecDt[i], te, plotConfig);
    }
    popLogStage();

    pushLogStage("CheckResult");
    auto errors = cutCellError(vecDisk, exactDisk, vecBox, plotConfig.range);
    printCellError(errors);

    t.~Timer();
    ::Timer::printStatistics();
    popLogStage();
    RecorderFinalize();
  }
}

TEST_CASE("Disk 5 deformation T = 2/4, Convergence Test.",
          "[Disk][Deformation][MARSn2D][Convergence]") {
  const auto* testName = "Disk5Deformation2";
  // const auto* testName = "Disk5Deformation4";
  ::Timer t(testName);
  auto dir = rootDir + "/results/TrackInterface/" + testName + "/";
  mkdir(dir.c_str(), 0755);
  RecorderInitialize(RecorderInfo{DebugLevel::INFO, dir});

  SECTION("4-th order disk") {
    constexpr int Order = 4;
    pushLogStage("Initialize");
    auto [vecDisk, exactDisk, radius, exactArea, exactLength, vecBox, vecN,
          aimOrder, vecHL, rTiny, nGrid, curvConfig, plotConfig, printDetail,
          t0, vecDt, te, timeIntegrator, velocityPtr] =
        diskTEST<Order>(rootDir + "/test/config/" + testName + ".json");
    popLogStage();

    pushLogStage("TrackInterface");
    for (uint i = 0; i < nGrid; ++i) {
      plotConfig.fName = to_string(Order) + "Circle" + "_grid" +
                         to_string(vecN[i]) + "_T" + to_string((int)te),
      plotConfig.fName = dir + plotConfig.fName;
      MARSn2D<Order, VectorFunction> CM(timeIntegrator, vecHL[i], rTiny,
                                        curvConfig, printDetail);

      CM.trackInterface(*velocityPtr, vecDisk[i], t0, vecDt[i], te, plotConfig);
    }
    popLogStage();

    pushLogStage("CheckResult");
    auto errors = cutCellError(vecDisk, exactDisk, vecBox, plotConfig.range);
    printCellError(errors);

    t.~Timer();
    ::Timer::printStatistics();
    popLogStage();
    RecorderFinalize();
  }
}