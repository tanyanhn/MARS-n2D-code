#pragma once

#include "InterfaceTracking/TimeIntegrator.h"
#include "Marsn2D/Elements.h"
#include "Marsn2D/InterfaceGraph.h"

namespace Marsn2D {

template <int Order, template <int> class VelocityField>
class MARSn2D {
  template <class T>
  using vector = std::vector<T>;
  using YS = YinSet<DIM, Order>;
  using IG = approxInterfaceGraph<Order>;
  using Edge = IG::Edge;

  using TIPtr = std::shared_ptr<TimeIntegrator<DIM, VelocityField>>;
  using Point = Vertex;

 public:
  struct CurvatureAdaptionConfig {
    bool used = false;
    Real rCMin;
    Interval<1> rhoCRange;
    std::function<Real(Real)> sigma;
  };
  enum { NONE = 0, NORMAL = 1, CUTCELL = 2 };
  struct PlotConfig {
    int output = NONE;
    bool plotInner;
    int opStride;
    std::string fName;
    Box<DIM> box;
    Interval<DIM> range;
  };

 public:
  MARSn2D() = delete;

  MARSn2D(TIPtr TI, Real hL, Real rTiny = 0.1,
          CurvatureAdaptionConfig curvConfig = {}, bool print = false)
      : TI_(TI),
        curvConfig_(curvConfig),
        hL_(hL),
        rTiny_(rTiny),
        printDetail(print) {}

  void trackInterface(const VelocityField<DIM> &v, IG &ig, Real StartTime,
                      Real dt, Real EndTime, const PlotConfig& plotConfig) const;

 private:
  void timeStep(const VelocityField<DIM> &v, IG &ig, Real tn, Real dt) const;
  void discreteFlowMap(const VectorFunction<2> &v, EdgeMark &marks, Real tn,
                       Real dt) const;

  void plot(const IG &ig, int step, const PlotConfig &plotConfig) const;

  void locateLongEdges(
      EdgeMark &marks, const vector<Real> &hL,
      vector<std::pair<unsigned int, unsigned int>> &indices2Num,
      Real efficientOfHL) const;
  void locateTinyEdges(EdgeMark &marks, const vector<Real> &hL,
                       vector<unsigned int> &indices, Real efficientOfHL) const;

  // insert {num} marks in i-th edge.
  void insertMarks(
      const VelocityField<2> &v, Real preT, Real dt,
      const vector<std::pair<unsigned int, unsigned int>> &indices2Num,
      EdgeMark &marks, Edge &preEdge, vector<Real> &curv) const;
  // delete i-th edge by remove i+1-th mark.
  // won't constitute delete mark.
  void removeMarks(const vector<unsigned int> &indices, EdgeMark &marks,
                   vector<Real> &curv) const;

  // update hL
  void updateHL(const vector<Real> &curv, const EdgeMark &marks,
                vector<Real> &hL) const;

  bool checkMarks(const EdgeMark &marks, const vector<Real>& hL) const;

 private:
  TIPtr TI_;

  CurvatureAdaptionConfig curvConfig_;
  Real hL_;
  Real rTiny_;

  bool printDetail;
  int removeTimesMax = 5;  
  int insertTimesMax = 5;
};
}  // namespace MARSn2D
