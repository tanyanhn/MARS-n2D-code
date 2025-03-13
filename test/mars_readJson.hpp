#pragma once
#include "testHeader.h"
#include <nlohmann/json.hpp>



// 参数结构体定义
struct SimulationParams {
  // Time parameters
  struct {
      Real te;
      Real t0;
      Real Cr;
      Real uM;
  } time;
  
  // grid parameters
  struct {
    int N0;
    Point lo;
    Point hi;
    int aimOrder;
    Real hLCoefficient;
    Real rTiny;
    int nGrid;
  } grid;

  // Curvature adaption parameters
  struct {
      bool curvUsed;
      Real rCMin;
      Real rhoMin;
      Real rhoMax;
  } curvature;

  // plot parameters
  struct {
      int plotOutput;
      bool plotInner;
      Real plotStride;
      bool printDetail;
  } plotConfig;

  // Domain disk parameters
  struct {
      Point center;
      rVec radius;
      std::vector<Real> parts;
  } domain;
};

// JSON解析适配器
inline void from_json(const nlohmann::json& j, Point& p) {
  p[0] = j[0].get<Real>();
  p[1] = j[1].get<Real>();
}

inline void from_json(const nlohmann::json& j, SimulationParams& params) {
  // 解析时间参数
  const auto& time = j.at("time");
  params.time.te = time.at("te");
  params.time.t0 = time.at("t0");
  params.time.Cr = time.at("Cr");
  params.time.uM = time.at("uM");

  // 解析网格参数
  const auto& grid = j.at("grid");
  params.grid.N0 = grid.at("N0");
  params.grid.lo = grid.at("lo");
  params.grid.hi = grid.at("hi");
  params.grid.aimOrder = grid.at("aimOrder");
  params.grid.hLCoefficient = grid.at("hLCoefficient");
  params.grid.rTiny = grid.at("rTiny");
  params.grid.nGrid = grid.at("nGrid");

  // 解析曲率参数
  const auto& curvature = j.at("curvature_adaption");
  params.curvature.curvUsed = curvature.at("curvUsed");
  params.curvature.rCMin = curvature.at("rCMin");
  params.curvature.rhoMin = curvature.at("rhoMin");
  params.curvature.rhoMax = curvature.at("rhoMax");

  const auto& plot = j.at("plot");
  params.plotConfig.plotOutput = plot.at("plotOutput");
  params.plotConfig.plotInner = plot.at("plotInner");
  params.plotConfig.plotStride = plot.at("plotStride");
  params.plotConfig.printDetail = plot.at("printDetail");

  // 解析域参数
  const auto& domain = j.at("domain_disk");
  params.domain.center = domain.at("center");
  params.domain.radius = domain.at("radius");
  params.domain.parts = domain.at("parts").get<std::vector<Real>>();
}

inline SimulationParams read_config(const std::string& filename) {
  std::ifstream ifs(filename);
  nlohmann::json config = nlohmann::json::parse(ifs);
  return config.get<SimulationParams>();
}