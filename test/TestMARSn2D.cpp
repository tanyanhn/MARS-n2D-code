#include <nlohmann/json.hpp>

#include "InterfaceTracking/ComputeError.h"
#include "Marsn2D/MARSn2D.h"
#include "Recorder/Timer.h"
#include "testHeader.h"

void trackInterfaceTest(const string& testName,
                        DebugLevel debugLevel = DebugLevel::OFF) {
  mkdir((rootDir + "/results").c_str(), 0755);
  mkdir((rootDir + "/results/TrackInterface").c_str(), 0755);
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
      auto localPlotConfig = plotConfig;
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
    std::string resultFileName = dir + "00Result.txt";
    printCellError(errors, resultFileName);
  }
  ::Timer::printStatistics();
  popLogStage();
  RecorderFinalize();
}

TEST_CASE("Disk 4 vortex T = 4, Order = 4", "[Disk4][Vortex][T4][Order4]") {
  trackInterfaceTest("Disk4VortexT4Order4");
}
TEST_CASE("Disk 4 vortex T = 4, Order = 6", "[Disk4][Vortex][T4][Order6]") {
  trackInterfaceTest("Disk4VortexT4Order6");
}
TEST_CASE("Disk 4 vortex T = 4, Order = 8", "[Disk4][Vortex][T4][Order8]") {
  trackInterfaceTest("Disk4VortexT4Order8");
}

TEST_CASE("Disk 4 vortex T = 8, Order = 4", "[Disk4][Vortex][T8][Order4]") {
  trackInterfaceTest("Disk4VortexT8Order4");
}
TEST_CASE("Disk 4 vortex T = 8, Order = 6", "[Disk4][Vortex][T8][Order6]") {
  trackInterfaceTest("Disk4VortexT8Order6");
}
TEST_CASE("Disk 4 vortex T = 8, Order = 8", "[Disk4][Vortex][T8][Order8]") {
  trackInterfaceTest("Disk4VortexT8Order8");
}

TEST_CASE("Disk 4 vortex T = 12, Order = 4", "[Disk4][Vortex][T12][Order4]") {
  trackInterfaceTest("Disk4VortexT12Order4");
}
TEST_CASE("Disk 4 vortex T = 12, Order = 6", "[Disk4][Vortex][T12][Order6]") {
  trackInterfaceTest("Disk4VortexT12Order6");
}
TEST_CASE("Disk 4 vortex T = 12, Order = 8", "[Disk4][Vortex][T12][Order8]") {
  trackInterfaceTest("Disk4VortexT12Order8");
}

TEST_CASE("Disk 4 vortex T = 16, Order = 4", "[Disk4][Vortex][T16][Order4]") {
  trackInterfaceTest("Disk4VortexT16Order4");
}
TEST_CASE("Disk 4 vortex T = 16, Order = 6", "[Disk4][Vortex][T16][Order6]") {
  trackInterfaceTest("Disk4VortexT16Order6");
}
TEST_CASE("Disk 4 vortex T = 16, Order = 8", "[Disk4][Vortex][T16][Order8]") {
  trackInterfaceTest("Disk4VortexT16Order8");
}

TEST_CASE("Disk 5 deformation T = 2, Order = 4",
          "[Disk5][Deformation][T2][Order4]") {
  trackInterfaceTest("Disk5DeformationT2Order4");
}
TEST_CASE("Disk 5 deformation T = 2, Order = 6",
          "[Disk5][Deformation][T2][Order6]") {
  trackInterfaceTest("Disk5DeformationT2Order6");
}
TEST_CASE("Disk 5 deformation T = 2, Order = 8",
          "[Disk5][Deformation][T2][Order8]") {
  trackInterfaceTest("Disk5DeformationT2Order8");
}

TEST_CASE("Disk 5 deformation T = 4, Order = 4",
          "[Disk5][Deformation][T4][Order4]") {
  trackInterfaceTest("Disk5DeformationT4Order4");
}
TEST_CASE("Disk 5 deformation T = 4, Order = 6",
          "[Disk5][Deformation][T4][Order6]") {
  trackInterfaceTest("Disk5DeformationT4Order6");
}
TEST_CASE("Disk 5 deformation T = 4, Order = 8",
          "[Disk5][Deformation][T4][Order8]") {
  trackInterfaceTest("Disk5DeformationT4Order8");
}

TEST_CASE("Disk 2 vortex T = 8, Order = 4", "[Disk2][Vortex][T8][Order4]") {
  trackInterfaceTest("Disk2VortexT8Order4");
}
TEST_CASE("Disk 2 Deformation T = 2, Order = 4",
          "[Disk2][Deformation][T2][Order4]") {
  trackInterfaceTest("Disk2DeformationT2Order4");
}

TEST_CASE("Disk 0 vortex T = 8, Order = 4", "[Disk0][Vortex][T8][Order4]") {
  trackInterfaceTest("Disk0VortexT8Order4");
}
TEST_CASE("Disk 0 vortex T = 12, Order = 4", "[Disk0][Vortex][T12][Order4]") {
  trackInterfaceTest("Disk0VortexT12Order4");
}

TEST_CASE("Graph 4.1, vortex T = 4, Order = 4", "[Graph41][Vortex][T4][Order4]") {
  constexpr int Order = 4;
  auto dir = rootDir + "/results/TrackInterface/Graph41VortexT4Order4/";
  mkdir(dir.c_str(), 0755);
  auto graph = InterfaceGraphFactory::createGraph41<Order>(0.001);

  auto yinsets = graph.approxYinSet();
  std::ofstream of(rootDir + "/results/TrackInterface/Graph41VortexT4Order4/"
                   "test.dat");
                   int num = yinsets.size();
  of.write((char*)&num, sizeof(int));
  for (int i = 0; i < yinsets.size(); ++i) yinsets[num - i - 1].dump(of);
}