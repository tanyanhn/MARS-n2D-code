#ifndef _MARS2DIMV_H_
#define _MARS2DIMV_H_

#include "MARS.h"
#include "Core/Curve.h"
#include "Core/Interval.h"
#include <vector>
#include <list>

template <int Order>
class MARS2DIMV : public MARS<2, Order>
{
    template <class T>
    using Vector = std::vector<T>;

    using YS = YinSet<2, Order>;

    using Point = Vec<Real, 2>;

    using Crv = Curve<2, Order>;

    using Base = MARS<2, Order>;

public:
    MARS2DIMV() = delete;

    MARS2DIMV(TimeIntegrator<2> *_TI, Real hL, Real rtiny = 0.1) : Base(_TI), chdLenRange(Vec<Real, 1>(rtiny * hL), Vec<Real, 1>(hL))
    {
    }

    void timeStep(const VectorFunction<2> &v, YS &ys, Real tn, Real dt);

private:
    void discreteFlowMap(const VectorFunction<2> &v, Vector<Point> &pts, Real tn, Real dt);

    Vector<unsigned int> removeSmallEdges(Vector<Point> &pts);

    Vector<unsigned int> splitLongEdges(const VectorFunction<2> &v, Vector<Point> &pts, const Crv &crv, Real tn, Real dt);

private:
    Interval<1> chdLenRange;
};

#endif