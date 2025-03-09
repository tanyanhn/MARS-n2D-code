#ifndef _MARS_H_
#define _MARS_H_

#include <string>

#include "Marsn2D/InterfaceGraph.h"
#include "TimeIntegrator.h"
#include "YinSet/YinSet.h"

template <int Dim, int Order, template <int> class VelocityField>
class MARS {
  using YS = YinSet<Dim, Order>;
  using IG = MARSn2D::approxInterfaceGraph<Order>;

 protected:
  TimeIntegrator<Dim, VelocityField> *TI;

 public:
  explicit MARS(TimeIntegrator<Dim, VelocityField> *TI_) : TI(TI_) {}

  virtual void timeStep(const VelocityField<Dim> &v, YS &ys, Real tn,
                        Real dt) = 0;

  virtual void trackInterface(const VelocityField<Dim> &v, YS &ys,
                              Real StartTime, Real dt, Real EndTime,
                              bool output, const std::string &fName,
                              int opStride);

  // Add for MARn2D.
  virtual void timeStep(const VelocityField<Dim> &v, IG &ig, Real tn, Real dt);

  virtual void trackInterface(const VelocityField<Dim> &v, IG &ig,
                              Real StartTime, Real dt, Real EndTime,
                              bool output, const std::string &fName,
                              int opStride);
};

#endif