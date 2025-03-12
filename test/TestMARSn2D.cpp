#include <nlohmann/json.hpp>

#include "InterfaceTracking/ERK.h"
#include "InterfaceTracking/RKButcher.h"
#include "InterfaceTracking/VelocityField.h"
#include "Marsn2D/MARSn2D.h"
#include "mars_readJson.hpp"
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
  // time
  const Real te = params.time.te;
  const Real t0 = params.time.t0;
  const Real Cr = params.time.Cr;
  const Real uM = params.time.uM;
  // grid
  const int N = params.grid.N;
  const rVec lo = params.grid.lo;
  const rVec hi = params.grid.hi;
  const int aimOrder = params.grid.aimOrder;
  const Real hLCoefficient = params.grid.hLCoefficient;
  const Real rTiny = params.grid.rTiny;
  const rVec h = (hi - lo) / N;
  const Real hL = hLCoefficient * std::pow(h[0] * h[1], 0.5 * aimOrder / Order);
  const Box<2> box(0, N - 1);
  const Real dt = std::min(h[0], h[1]) / uM * Cr;
  const auto timeIntegrator =
      make_shared<ERK<2, RK::ClassicRK4, VectorFunction>>();
  const Real T = te;
  const Vortex vortex(T);

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
    .fName = to_string(Order) + "Circle" + "_grid" + to_string(N) + "_T" +
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
  const auto disk =
      Generator::createDiskGraph<Order>(center, radius, hL, parts);

  return make_tuple(disk, radius, exactArea, exactLength, box, N, h, hL, rTiny,
                    curvConfig, plotConfig, printDetail, t0, dt, te, timeIntegrator, vortex);
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

TEST_CASE("Disk 4 vortex4, RK4.", "[Disk][Vortex][4][MARSn2D]") {
  const auto* testName = "Disk4Vortex4";
  auto dir = rootDir + "/results/TrackInterface/" + testName + "/";
  mkdir(dir.c_str(), 0755);
  SECTION("2-th order disk") {
    constexpr int Order = 2;
    auto [disk, radius, exactArea, exactLength, box, N, h, hL, rTiny,
          curvConfig, plotConfig, printDetail, t0, dt, te, timeIntegrator,
          vortex] =
        diskTEST<Order>(rootDir + "/test/config/" + testName + ".json");
    plotConfig.fName = dir + plotConfig.fName;
    MARSn2D<Order, VectorFunction> CM(timeIntegrator, hL, rTiny, curvConfig,
      printDetail);

    CM.trackInterface(vortex, disk, t0, dt, te, plotConfig);

    INFO("Integral Area and length.");
    auto yinSets = disk.approxYinSet();

    checkResult<Order>(yinSets, box, plotConfig.range, addInner, radius, N, h,
                       output, dir, plotConfig.fName, exactLength, exactArea);
  }
}