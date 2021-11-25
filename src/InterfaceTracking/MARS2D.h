#ifndef _MARS2D_
#define _MARS2D_

#include "MARS.h"
#include "Core/Vec.h"
#include "Core/Curve.h"
#include "Core/Interval.h"
#include <vector>


//Real tol = 1e-15;

template <class T>
using Vector = std::vector<T>;

template <int Order, template <typename...> class Container>
class MARS2D : public MARS<2, Order>
{
    using YS = YinSet<2, Order>;

    using Point = Vec<Real, 2>;

    using Crv = Curve<2, Order>;

    using Base = MARS<2, Order>;

private:
    Interval<1> chdLenRange;

    void discreteFlowMap(const IT_VectorFunction<2> &v, Container<Point> &pts, Real tn, Real dt);

    Vector<unsigned int> removeSmallEdges(Container<Point> &pts);

    Vector<unsigned int> splitLongEdges(const IT_VectorFunction<2> &v, Container<Point> &pts, const Crv &crv, Real tn, Real dt);

public:
    MARS2D() = delete;
    MARS2D(IT_TimeIntegrator<2> *_TI, Real hL, Real rtiny = 0.1) : Base(_TI), chdLenRange(Vec<Real, 1>(rtiny * hL), Vec<Real, 1>(hL))
    {
    }

    void timeStep(const IT_VectorFunction<2> &v, YS &ys, Real tn, Real dt);
};

#endif