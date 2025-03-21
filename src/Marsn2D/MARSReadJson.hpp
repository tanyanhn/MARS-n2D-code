#pragma once
// #include "testHeader.h"
#include <nlohmann/json.hpp>

#include "Core/Config.h"
#include "Core/Vec.h"
#include "InterfaceTracking/VelocityField.h"

using Point = Vec<Real, 2>;
using rVec = Vec<Real, 2>;

// 参数结构体定义
struct SimulationParams {
  // Tolerance parameters
  struct {
    Real distTol;
    Real newtonTol;
    int newtonMaxIter;
  } tolerance;
  // Time parameters
  struct {
    Real te;
    Real t0;
    Real Cr;
    Real uM;
  } time;

  // field parameters
  struct {
    std::shared_ptr<VectorFunction<2>> velocity;
    std::string name;
    ParamWrapper pars;
  } field;

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
    int plotStride;
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
  const auto& tolerance = j.at("tolerance");
  params.tolerance.distTol = tolerance.at("distTol").get<Real>();
  params.tolerance.newtonTol = tolerance.at("newtonTol").get<Real>();
  params.tolerance.newtonMaxIter = tolerance.at("newtonMaxIter").get<int>();
  // 解析时间参数
  const auto& time = j.at("time");
  params.time.te = time.at("te");
  params.time.t0 = time.at("t0");
  params.time.Cr = time.at("Cr");
  params.time.uM = time.at("uM");

  const auto& field = j.at("field");
  params.field.name = field.at("velocityName");
  auto constructPars = [](Real te, const nlohmann::json& j) {
    if (j.empty()) return ParamWrapper(te);
    if (j.is_array() && j[0].is_number_integer())
      return ParamWrapper(te, j[0].get<int>()); // Deformation
    if (j.is_array() && j[0].is_number_float())
      return ParamWrapper(te, j[0].get<Real>());  // Vortex
    return ParamWrapper(te); // unknow default;
  };
  params.field.pars = constructPars(params.time.te, field.at("pars"));
  params.field.velocity.reset(Vector2DFactory::getInstance().create(
      params.field.name, params.field.pars));

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
  if (!ifs.is_open()) {
    std::cerr << "Error: Failed to open config file: " << filename << std::endl;
    exit(1);
  }
  nlohmann::json config = nlohmann::json::parse(ifs);
  return config.get<SimulationParams>();
}