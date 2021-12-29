#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "Core/Config.h"
#include "Core/Vec.h"
#include "VectorFunction.h"
#include <vector>

template <int Dim>
class TimeIntegrator
{
    template <class T>
    using Vector = std::vector<T>;

    using Point = Vec<Real, Dim>;

public:
    virtual ~TimeIntegrator() {}

    virtual void timeStep(const VectorFunction<Dim> &v, Point &pt, Real tn, Real k) = 0;

    virtual void timeStep(const VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real k) = 0;
};

#endif