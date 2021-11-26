#ifndef _MARS2D_H_
#define _MARS2D_H_

#include "MARS.h"
#include "Core/Curve.h"
#include "Core/Interval.h"
#include <vector>

template <int Order, template <typename...> class Container>
class MARS2D : public MARS<2, Order>
{
    template <class T>
    using Vector = std::vector<T>;

    using YS = YinSet<2, Order>;

    using Point = Vec<Real, 2>;

    using Crv = Curve<2, Order>;

    using Base = MARS<2, Order>;

public:
    MARS2D() = delete;

    MARS2D(TimeIntegrator<2> *_TI, Real hL, Real rtiny = 0.1) : Base(_TI), chdLenRange(Vec<Real, 1>(rtiny * hL), Vec<Real, 1>(hL))
    {
    }

    void timeStep(const VectorFunction<2> &v, YS &ys, Real tn, Real dt);

private:
    void discreteFlowMap(const VectorFunction<2> &v, Container<Point> &pts, Real tn, Real dt);

    Vector<unsigned int> removeSmallEdges(Container<Point> &pts);

    Vector<unsigned int> splitLongEdges(const VectorFunction<2> &v, Container<Point> &pts, const Crv &crv, Real tn, Real dt);

private:
    Interval<1> chdLenRange;
};

#endif