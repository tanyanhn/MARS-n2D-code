#ifndef _MARS_H_
#define _MARS_H_

#include <string>

#include "TimeIntegrator.h"
#include "YinSet/YinSet.h"

template <int Dim, int Order, template <int> class VelocityField>
class MARS {
  using YS = YinSet<Dim, Order>;

 protected:
  TimeIntegrator<Dim, VelocityField> *TI;

 public:
  MARS(TimeIntegrator<Dim, VelocityField> *_TI) : TI(_TI) {}

  virtual void timeStep(const VelocityField<Dim> &v, YS &ys, Real tn,
                        Real dt) = 0;

  virtual void trackInterface(const VelocityField<Dim> &v, YS &ys,
                              Real StartTime, Real dt, Real EndTime,
                              bool output = false, const std::string& fName = "",
                              int opstride = 20);
};

#endif