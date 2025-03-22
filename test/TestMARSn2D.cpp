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
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
auto diskTEST(const std::string& jsonFile) {
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
  // const rVec h0 = (hi - lo) / N0;
  // const Real hL0 =
  //     hLCoefficient * std::pow(h0[0] * h0[1], 0.5 * aimOrder / Order);
  // const Box<DIM> box(0, N0 - 1);
  // const Real dt0 = std::min(h0[0], h0[1]) / uM * Cr;
  std::shared_ptr<TimeIntegrator<DIM, VectorFunction>> timeIntegrator;
  if (aimOrder == 4) {
    timeIntegrator = make_shared<ERK<DIM, RK::ClassicRK4, VectorFunction>>();
  } else if (aimOrder == 6) {
    timeIntegrator = make_shared<ERK<DIM, RK::Verner6, VectorFunction>>();
  } else if (aimOrder == 8) {
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
  const int plotStride = params.plotConfig.plotStride;
  const bool printDetail = params.plotConfig.printDetail;
  const struct MARSn2D<Order, VectorFunction>::PlotConfig plotConfig {
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
        Generator::createDiskGraph<Order>(center, radius, hL / 2, parts);
    vecDisk.emplace_back(std::move(disk));
    vecN.emplace_back(N);
    vecBox.emplace_back(0, N - 1);
    // vecH.emplace_back(h);
    vecHL.emplace_back(hL);
    vecDt.emplace_back(dt);
    N = N * 2;
  }

  // accurate solution
  Real hL = vecHL.back() / 10;
  auto exactDisk =
      Generator::createDiskGraph<Order>(center, radius, hL / 2, parts);

  return make_tuple(vecDisk, exactDisk, radius, exactArea, exactLength, vecBox,
                    vecN, aimOrder, vecHL, rTiny, nGrid, curvConfig, plotConfig,
                    printDetail, t0, vecDt, te, timeIntegrator, velocityPtr);
  // }
}

void trackInterfaceTest(const string& testName,
                        DebugLevel debugLevel = DebugLevel::OFF) {
  auto dir = rootDir + "/results/TrackInterface/" + testName + "/";
  mkdir(dir.c_str(), 0755);
  RecorderInitialize(RecorderInfo{debugLevel, dir});
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
      // 每个线程创建自己的plotConfig副本
      auto localPlotConfig = plotConfig;  // 假设PlotConfig可复制
      localPlotConfig.fName =
          to_string(Order) + "Circle" + "_grid" + to_string(vecN[i]);
      localPlotConfig.fName = dir + localPlotConfig.fName;
      MARSn2D<Order, VectorFunction> CM(timeIntegrator, vecHL[i], rTiny,
                                        curvConfig, printDetail);
      CM.trackInterface(*velocityPtr, vecDisk[i], t0, vecDt[i], te,
                        localPlotConfig);
    }
    popLogStage();
    pushLogStage("CheckResult");
    auto errors = cutCellError(vecDisk, exactDisk, vecBox, plotConfig.range);
    printCellError(errors);
  }
  ::Timer::printStatistics();
  popLogStage();
  RecorderFinalize();
}

TEST_CASE("Disk 4 vortex T = 4, Order = 4", "[Disk][Vortex][T4][Order4]") {
  trackInterfaceTest("Disk4VortexT4Order4");
}
TEST_CASE("Disk 4 vortex T = 4, Order = 6", "[Disk][Vortex][T4][Order6]") {
  trackInterfaceTest("Disk4VortexT4Order6");
}
TEST_CASE("Disk 4 vortex T = 4, Order = 8", "[Disk][Vortex][T4][Order8]") {
  trackInterfaceTest("Disk4VortexT4Order8");
}

TEST_CASE("Disk 4 vortex T = 8, Order = 4", "[Disk][Vortex][T8][Order4]") {
  trackInterfaceTest("Disk4VortexT8Order4");
}
TEST_CASE("Disk 4 vortex T = 8, Order = 6", "[Disk][Vortex][T8][Order6]") {
  trackInterfaceTest("Disk4VortexT8Order6");
}
TEST_CASE("Disk 4 vortex T = 8, Order = 8", "[Disk][Vortex][T8][Order8]") {
  trackInterfaceTest("Disk4VortexT8Order8");
}

TEST_CASE("Disk 4 vortex T = 12, Order = 4", "[Disk][Vortex][T12][Order4]") {
  trackInterfaceTest("Disk4VortexT12Order4");
}
TEST_CASE("Disk 4 vortex T = 12, Order = 6", "[Disk][Vortex][T12][Order6]") {
  trackInterfaceTest("Disk4VortexT12Order6");
}
TEST_CASE("Disk 4 vortex T = 12, Order = 8", "[Disk][Vortex][T12][Order8]") {
  trackInterfaceTest("Disk4VortexT12Order8");
}

TEST_CASE("Disk 4 vortex T = 16, Order = 4", "[Disk][Vortex][T16][Order4]") {
  trackInterfaceTest("Disk4VortexT16Order4");
}
TEST_CASE("Disk 4 vortex T = 16, Order = 6", "[Disk][Vortex][T16][Order6]") {
  trackInterfaceTest("Disk4VortexT16Order6");
}
TEST_CASE("Disk 4 vortex T = 16, Order = 8", "[Disk][Vortex][T16][Order8]") {
  trackInterfaceTest("Disk4VortexT16Order8");
}

TEST_CASE("Disk 5 deformation T = 2, Order = 4",
          "[Disk][Deformation][T2][Order4]") {
  trackInterfaceTest("Disk5DeformationT2Order4");
}
TEST_CASE("Disk 5 deformation T = 2, Order = 6",
          "[Disk][Deformation][T2][Order6]") {
  trackInterfaceTest("Disk5DeformationT2Order6");
}
TEST_CASE("Disk 5 deformation T = 2, Order = 8",
          "[Disk][Deformation][T2][Order8]") {
  trackInterfaceTest("Disk5DeformationT2Order8");
}

TEST_CASE("Disk 5 deformation T = 4, Order = 4",
          "[Disk][Deformation][T4][Order4]") {
  trackInterfaceTest("Disk5DeformationT4Order4");
}
TEST_CASE("Disk 5 deformation T = 4, Order = 6",
          "[Disk][Deformation][T4][Order6]") {
  trackInterfaceTest("Disk5DeformationT4Order6");
}
TEST_CASE("Disk 5 deformation T = 4, Order = 8",
          "[Disk][Deformation][T4][Order8]") {
  trackInterfaceTest("Disk5DeformationT4Order8");
}
