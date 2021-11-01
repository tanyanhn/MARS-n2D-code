#ifndef _CUBICMARS_
#define _CUBICMARS_

#include "InterfaceTracking/TimeIntegrator.h"
#include "InterfaceTracking/ExplicitRK.h"
#include "YinSet/YinSet.h"
#include "YinSet/SegmentedRealizableSpadjor.h"
#include "Core/Vec.h"
#include "Core/Curve.h"
#include "Core/Polynomial.h"
#include "Core/Interval.h"
#include <cmath>
#include <vector>
#include <functional>

Real tol = 1e-15;

template <class T>
using Vector = std::vector<T>;

template <int Dim, int Order>
class MARS
{
    using YS = YinSet<Dim, Order>;

    using Point = Vec<Real, Dim>;

    using Crv = Curve<Dim, Order>;

private:
    IT_TimeIntegrator<Dim> *TI;
    Interval<1> chdLenRange;

public:
    MARS() = delete;
    MARS(IT_TimeIntegrator<Dim> *_TI, Real hL, Real rtiny = 0.1) : TI(_TI), chdLenRange(Vec<Real, 1>(rtiny * hL), Vec<Real, 1>(hL))
    {
    }

    void discreteFlowMap(const IT_VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real dt);

    Vector<unsigned int> removeSmallEdges(Vector<Point> &pts);

    Vector<unsigned int> splitLongEdges(const IT_VectorFunction<Dim> &v, Vector<Point> &pts, Crv &crv, Real tn, Real dt);

    void timeStep(const IT_VectorFunction<Dim> &v, YS &ys, Real tn, Real dt);

    void trackInterface(const IT_VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime);
};

template <int Dim, int Order>
void MARS<Dim, Order>::discreteFlowMap(const IT_VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real dt)
{
    TI->timeStep(v, pts, tn, dt);
    return;
}

template <int Dim>
bool remove(Vector<unsigned int> &ids, const Vector<Vec<Real, Dim>> &pts, Real lowBound)
{
    size_t num = ids.size();
    Vector<Real> dist(num - 1);
    //dist.resize(num - 1);
    for (size_t i = 0; i < num - 1; i++)
    {
        dist[i] = norm(pts[ids[i + 1]] - pts[ids[i]]);
    }

    //mark weather the pre-node is candidate to be erased.
    bool predelete = false;
    auto it = ids.begin();
    ++it;

    size_t i = 1;
    if (dist[0] < lowBound)
    {
        predelete = true;
        it = ids.erase(it);
        i++;
    }

    while (i < num - 2)
    {
        if (dist[i] >= lowBound)
        {
            predelete = false;
            ++it;
            i++;
        }
        else
        {
            if (dist[i - 1] <= dist[i + 1])
            {
                if (predelete == false)
                {
                    it = ids.erase(it);
                    i++;
                    predelete = true;
                }
                else
                {
                    ++it;
                    i++;
                    predelete = false;
                }
            }
            else
            {
                it = ids.erase(it + 1);
                i += 2;
                predelete = true;
            }
        }
    }
    if (i == num - 1)
    {
        return num != ids.size();
    }
    else
    {
        if (predelete == true || dist[i] >= lowBound)
            return num != ids.size();
        it = ids.erase(it);
        return true;
    }
}

template <int Dim, int Order>
Vector<unsigned int> MARS<Dim, Order>::removeSmallEdges(Vector<Point> &pts)
{
    size_t Num = pts.size();
    Vector<unsigned int> ids(Num);
    for (size_t i = 0; i < Num; i++)
    {
        ids[i] = i;
    }

    while (true)
    {
        if (!remove(ids, pts, (chdLenRange.lo())[0]))
            break;
    }
    Vector<Point> npts;
    npts.resize(ids.size());
    for (size_t i = 0; i < ids.size(); i++)
    {
        npts[i] = pts[ids[i]];
    }
    pts = npts;
    return ids;
}

template <int Dim, int Order>
Vector<unsigned int> MARS<Dim, Order>::splitLongEdges(const IT_VectorFunction<Dim> &v, Vector<Point> &pts, Crv &crv, Real tn, Real dt)
{
    assert(crv.isClosed(tol));
    size_t Num = pts.size();
    Vector<Point> oldpts(Num);
    Vector<unsigned int> ids;
    Vector<Polynomial<Order, Point>> poly = crv.getPolys();
    Vector<Real> knots = crv.getKnots();

    for (size_t i = 0; i < Num - 1; i++)
    {
        oldpts[i] = (poly[i])(0);
    }
    oldpts[Num - 1] = (poly[0])(0);

    std::function<void(Real, Real, Vector<Real> &)> split;

    split = [&](Real left, Real right, Vector<Real> &chdt) -> void
    {
        Point lpt = crv(left);
        lpt = TI->timeStep(v, lpt, tn, dt);
        Point rpt = crv(right);
        rpt = TI->timeStep(v, rpt, tn, dt);
        if (norm(rpt - lpt) > (chdLenRange.hi())[0])
        {
            int N = (int)ceil(norm(rpt - lpt) / (chdLenRange.hi())[0]);
            Real dl = (right - left) / N;
            for (int j = 0; j < N - 1; j++)
            {
                split(left + dl * j, left + dl * (j + 1), chdt);
                chdt.push_back(left + dl * (j + 1));
            }
            split(right - dl, right, chdt);
        }
        else
        {
            return;
        }
    };

    auto it = pts.begin();
    auto oit = oldpts.begin();
    ++it;
    ++oit;
    int count = 1;

    Real dist;
    for (size_t i = 0; i < Num - 1; i++)
    {
        dist = norm(pts[i + 1] - pts[i]);
        if (dist <= (chdLenRange.hi())[0])
        {
            ++it;
            ++oit;
            count++;
            continue;
        }
        Vector<Real> chdt;
        split(knots[i], knots[i + 1], chdt);
        size_t sz = ids.size();
        ids.resize(sz + chdt.size());
        for (size_t j = sz; j < ids.size(); j++)
        {
            count++;
            ids[j] = count;
        }
        count++;
        for (size_t j = chdt.size(); j >= 1; j--)
        {
            Point opt = crv(chdt[j-1]);
            oit = oldpts.emplace(oit, opt);
            it = pts.emplace(it, TI->timeStep(v, opt, tn, dt));
        }
        it += chdt.size() + 1;
        oit += chdt.size() + 1;
    }
    crv = fitCurve<Order>(oldpts, true);
    return ids;
}

template <int Dim, int Order>
void MARS<Dim, Order>::timeStep(const IT_VectorFunction<Dim> &v, YS &ys, Real tn, Real dt)
{
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    for (auto &crv : vcrv)
    {
        assert(crv.isClosed(tol));

        //pts: next time's points; oldpts: now time's points
        Vector<Point> oldpts;
        Vector<Polynomial<Order,Point>> poly = crv.getPolys();
        oldpts.resize(poly.size() + 1);
        for (size_t i = 0; i < poly.size(); i++)
        {
            oldpts[i] = (poly[i])[0];
        }
        oldpts[poly.size()] = (poly[0])[0];

        //get the next time's points
        Vector<Point> pts = oldpts;
        TI->timeStep(v, pts, tn, dt);

        //split long edges, and regenerate the old curve
        Vector<unsigned int> splitids = splitLongEdges(v, pts, crv, tn, dt);

        //remove small edges, update pts and old pts simultaneously
        Vector<unsigned int> removeids = removeSmallEdges(pts);

        crv = fitCurve<Order>(pts, true);
    }
    ys = YS(SegmentedRealizableSpadjor<Order>(vcrv), tol);
    return;
}

template <int Dim, int Order>
void MARS<Dim, Order>::trackInterface(const IT_VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime)
{
    Real T = StartTime;
    Real t = dt;
    while (T < EndTime)
    {
        if (EndTime - T < dt)
        {
            t = EndTime - T;
        }
        timeStep(v, ys, T, t);
        T += t;
    }
    return;
}

template class MARS<2, 4>;

#endif