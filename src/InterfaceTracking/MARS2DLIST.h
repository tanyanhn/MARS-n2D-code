#ifndef _MARS2DLIST_H_
#define _MARS2DLIST_H_

#include "MARS.h"
#include "Core/Curve.h"
#include "Core/Interval.h"
#include <vector>
#include <list>

template <int Order>
class MARS2DLIST : public MARS<2, Order>
{
    template <class T>
    using List = std::list<T>;

    template <class T>
    using Vector = std::vector<T>;

    using YS = YinSet<2, Order>;

    using Point = Vec<Real, 2>;

    using Crv = Curve<2, Order>;

    using Base = MARS<2, Order>;

public:
    MARS2DLIST() = delete;

    MARS2DLIST(TimeIntegrator<2> *_TI, Real hL, Real rtiny = 0.1) : Base(_TI), chdLenRange(Vec<Real, 1>(rtiny * hL), Vec<Real, 1>(hL))
    {
    }

    void timeStep(const VectorFunction<2> &v, YS &ys, Real tn, Real dt){}

    void timeStep(const VectorFunction<2> &v, YS &ys, Vector<List<Point>> &vpts, Real tn, Real dt);

    void trackInterface(const VectorFunction<2> &v, YS &ys, Real StartTime, Real dt, Real EndTime, bool output = false, std::string fName = "", int opstride = 20);

private:
    void discreteFlowMap(const VectorFunction<2> &v, List<Point> &pts, Real tn, Real dt);

    int removeSmallEdges(List<Point> &pts);

    int splitLongEdges(const VectorFunction<2> &v, List<Point> &pts, const Crv &crv, Real tn, Real dt);

private:
    Interval<1> chdLenRange;
};

#endif