#ifndef _MARS_H_
#define _MARS_H_

#include "TimeIntegrator.h"
#include "YinSet/YinSet.h"
#include <string>

template <int Dim, int Order>
class MARS
{
    using YS = YinSet<Dim, Order>;

protected:
    IT_TimeIntegrator<Dim> *TI;

public:
    MARS(IT_TimeIntegrator<Dim> *_TI) : TI(_TI) {}

    virtual void timeStep(const IT_VectorFunction<Dim> &v, YS &ys, Real tn, Real dt) = 0;

    virtual void trackInterface(const IT_VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime, bool output = false, std::string fName = "", int opstride = 20);
};

#endif