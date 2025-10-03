#include <nlohmann/json.hpp>

#include "InterfaceTracking/ComputeError.h"
#include "Marsn2D/MARSn2D.h"
#include "Recorder/Timer.h"
#include "testHeader.h"

OPTNONE_FUNC
void trackInterfaceTest(const string& testName,
                        DebugLevel debugLevel = DebugLevel::OFF) {
  mkdir((rootDir + "/results").c_str(), 0755);
  mkdir((rootDir + "/results/TrackInterface").c_str(), 0755);
  auto dir = rootDir + "/results/TrackInterface/" + testName + "/";
  mkdir(dir.c_str(), 0755);
  RecorderInitialize(RecorderInfo{debugLevel, dir});
  // ::Timer::initial();
  SECTION("4-th order") {
    constexpr int Order = 4;
    pushLogStage("Initialize");
    auto [vecDomain, exactDomain, radius, exactArea, exactLength, vecBox, vecN,
          aimOrder, vecHL, rTiny, nGrid, curvConfig, plotConfig, printDetail,
          t0, vecDt, te, T, timeIntegrator, velocityPtr] =
        diskTEST<Order>(rootDir + "/test/config/" + testName + ".json");
    popLogStage();
    pushLogStage("TrackInterface");
    for (uint i = 0; i < nGrid; ++i) {
      auto localPlotConfig = plotConfig;
      localPlotConfig.fName = std::string("1Order") + to_string(Order) +
                              "_grid" + to_string(vecN[i]);
      localPlotConfig.fName = dir + localPlotConfig.fName;
      localPlotConfig.box = vecBox[i];
      MARSn2D<Order, VectorFunction> CM(timeIntegrator, vecHL[i], rTiny,
                                        curvConfig, printDetail);
      CM.trackInterface(*velocityPtr, vecDomain[i], t0, vecDt[i], te,
                        localPlotConfig);
    }
    popLogStage();
    pushLogStage("CheckResult");
    // std::ofstream of(dir + "00test.dat");
    // for (auto& domain : vecDomain) {
    //   auto yinSets = domain.approxYinSet();
    //   int num = yinSets.size();
    //   of.write((char*)&num, sizeof(int));
    //   for (auto& yin : yinSets) yin.dump(of);
    //   of.close();
    // }
    auto errors =
        cutCellError(vecDomain, exactDomain, vecBox, plotConfig.range, te != T);
    std::string resultFileName = dir + "00Result.txt";
    printCellError(errors, resultFileName);
  }
  ::Timer::printStatistics();
  popLogStage();
  RecorderFinalize();
}

TEST_CASE("Dynamic Graph 4.1, vortex T = 16, Order = 4",
          "[Graph41][Vortex][T16][Order4][Dynamic]") {
  trackInterfaceTest("Graph41VortexT16Order4Dynamic");
}

TEST_CASE("Dynamic Graph 4.1, Deformation T = 4, Order = 4",
          "[Graph41][Deformation][T4][Order4][Dynamic]") {
  trackInterfaceTest("Graph41DeformationT4Order4Dynamic");
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
          "[Disk5][Deformation][T2][Order4][Compare][origin]") {
  trackInterfaceTest("Disk5DeformationT2Order4");
}
TEST_CASE("Disk 5 deformation T = 2, Order = 4, Compare hL",
          "[Disk5][Deformation][T2][Order4][Compare][hL]") {
  trackInterfaceTest("Disk5DeformationT2Order4HL");
}
TEST_CASE("Disk 5 deformation T = 2, Order = 4, Compare k = h",
          "[Disk5][Deformation][T2][Order4][Compare][k=h]") {
  trackInterfaceTest("Disk5DeformationT2Order4K=h");
}
TEST_CASE("Disk 5 deformation T = 2, Order = 4, Compare Richardson",
          "[Disk5][Deformation][T2][Order4][Compare][Richardson]") {
  trackInterfaceTest("Disk5DeformationT2Order4Richard");
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

TEST_CASE("Graph 4.1, vortex T = 4, Order = 4",
          "[Graph41][Vortex][T4][Order4]") {
  trackInterfaceTest("Graph41VortexT4Order4");
}
TEST_CASE("Graph 4.1, vortex T = 4, Order = 6",
          "[Graph41][Vortex][T4][Order6]") {
  trackInterfaceTest("Graph41VortexT4Order6");
}
TEST_CASE("Graph 4.1, vortex T = 4, Order = 8",
          "[Graph41][Vortex][T4][Order8]") {
  trackInterfaceTest("Graph41VortexT4Order8");
}

TEST_CASE("Graph 4.1, vortex T = 8, Order = 4",
          "[Graph41][Vortex][T8][Order4]") {
  trackInterfaceTest("Graph41VortexT8Order4");
}
TEST_CASE("Graph 4.1, vortex T = 8, Order = 6",
          "[Graph41][Vortex][T8][Order6]") {
  trackInterfaceTest("Graph41VortexT8Order6");
}
TEST_CASE("Graph 4.1, vortex T = 8, Order = 8",
          "[Graph41][Vortex][T8][Order8]") {
  trackInterfaceTest("Graph41VortexT8Order8");
}
TEST_CASE("Graph 4.1, vortex T = 12, Order = 4",
          "[Graph41][Vortex][T12][Order4]") {
  trackInterfaceTest("Graph41VortexT12Order4");
}
TEST_CASE("Graph 4.1, vortex T = 12, Order = 6",
          "[Graph41][Vortex][T12][Order6]") {
  trackInterfaceTest("Graph41VortexT12Order6");
}
TEST_CASE("Graph 4.1, vortex T = 12, Order = 8",
          "[Graph41][Vortex][T12][Order8]") {
  trackInterfaceTest("Graph41VortexT12Order8");
}
TEST_CASE("Graph 4.1, vortex T = 16, Order = 4",
          "[Graph41][Vortex][T16][Order4]") {
  trackInterfaceTest("Graph41VortexT16Order4");
}
TEST_CASE("Graph 4.1, vortex T = 16, Order = 6",
          "[Graph41][Vortex][T16][Order6]") {
  trackInterfaceTest("Graph41VortexT16Order6");
}
TEST_CASE("Graph 4.1, vortex T = 16, Order = 8",
          "[Graph41][Vortex][T16][Order8]") {
  trackInterfaceTest("Graph41VortexT16Order8");
}

TEST_CASE("Graph 4.1, Deformation T = 2, Order = 4",
          "[Graph41][Deformation][T2][Order4]") {
  trackInterfaceTest("Graph41DeformationT2Order4");
}
TEST_CASE("Graph 4.1, Deformation T = 2, Order = 6",
          "[Graph41][Deformation][T2][Order6]") {
  trackInterfaceTest("Graph41DeformationT2Order6");
}
TEST_CASE("Graph 4.1, Deformation T = 4, Order = 8",
          "[Graph41][Deformation][T2][Order8]") {
  trackInterfaceTest("Graph41DeformationT2Order8");
}
TEST_CASE("Graph 4.1, Deformation T = 4, Order = 4",
          "[Graph41][Deformation][T4][Order4]") {
  trackInterfaceTest("Graph41DeformationT4Order4");
}
TEST_CASE("Graph 4.1, Deformation T = 4, Order = 6",
          "[Graph41][Deformation][T4][Order6]") {
  trackInterfaceTest("Graph41DeformationT4Order6");
}
TEST_CASE("Graph 4.1, Deformation T = 4, Order = 8",
          "[Graph41][Deformation][T4][Order8]") {
  trackInterfaceTest("Graph41DeformationT4Order8");
}