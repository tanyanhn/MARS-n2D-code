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
#include <iomanip>
#include <algorithm>
#include <string>

//Real tol = 1e-15;

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

    Vector<unsigned int> splitLongEdges(const IT_VectorFunction<Dim> &v, Vector<Point> &pts, const Crv &crv, Real tn, Real dt);

    void timeStep(const IT_VectorFunction<Dim> &v, YS &ys, Real tn, Real dt);

    void trackInterface(const IT_VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime, bool output = false, std::string fName = "", int opstride = 20);
};
/*
template <int Dim, int Order>
void MARS<Dim, Order>::discreteFlowMap(const IT_VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real dt)
{
    TI->timeStep(v, pts, tn, dt);
    return;
}

template <int Dim>
bool remove(Vector<unsigned int> &ids, Vector<Vec<Real, Dim>> &pts, Real lowBound)
{
    int num = ids.size();
    Vector<Real> dist(num - 1);
    //dist.resize(num - 1);
    for (int i = 0; i < num - 1; i++)
    {
        dist[i] = norm(pts[i + 1] - pts[i]);
    }

    //mark weather the pre-node is candidate to be erased.
    bool predelete = false;
    auto it = ids.begin();
    auto pit = pts.begin();
    ++it;
    ++pit;

    int i = 1;
    if (dist[0] < lowBound)
    {
        predelete = true;
        it = ids.erase(it);
        pit = pts.erase(pit);
        i++;
    }

    while (i < num - 2)
    {
        if (dist[i] >= lowBound)
        {
            predelete = false;
            ++it;
            ++pit;
            i++;
        }
        else
        {
            if (dist[i - 1] <= dist[i + 1])
            {
                if (predelete == false)
                {
                    it = ids.erase(it);
                    pit = pts.erase(pit);
                    i++;
                    predelete = true;
                }
                else
                {
                    ++it;
                    ++pit;
                    i++;
                    predelete = false;
                }
            }
            else
            {
                it = ids.erase(it + 1);
                pit = pts.erase(pit + 1);
                i += 2;
                predelete = true;
            }
        }
    }
    if (i == num - 1)
    {
        return num != (int)ids.size();
    }
    else
    {
        if (predelete == true || dist[i] >= lowBound)
            return num != (int)ids.size();
        it = ids.erase(it);
        pit = pts.erase(pit);
        return true;
    }
}

template <int Dim, int Order>
Vector<unsigned int> MARS<Dim, Order>::removeSmallEdges(Vector<Point> &pts)
{
    int num = pts.size();
    Vector<unsigned int> ids(num);
    for (int i = 0; i < num; i++)
    {
        ids[i] = i;
    }

    while (true)
    {
        if (!remove(ids, pts, (chdLenRange.lo())[0]))
            break;
    }
    Vector<unsigned int> rids(num - ids.size());
    int idsid = 0;
    unsigned int allid = 0;
    int rid = 0;
    while (rid < (int)rids.size() - 1)
    {
        if (ids[idsid] != allid)
        {
            rids[rid] = allid;
            rid++;
            allid++;
        }
        else
        {
            idsid++;
            allid++;
        }
    }
    return rids;
}

template <int Dim, int Order>
Vector<unsigned int> MARS<Dim, Order>::splitLongEdges(const IT_VectorFunction<Dim> &v, Vector<Point> &pts, const Crv &crv, Real tn, Real dt)
{
    assert(crv.isClosed(tol));
    int num = pts.size();
    Vector<Point> oldpts(num);
    Vector<unsigned int> ids;
    Vector<Polynomial<Order, Point>> poly = crv.getPolys();
    Vector<Real> knots = crv.getKnots();

    for (int i = 0; i < num - 1; i++)
    {
        oldpts[i] = (poly[i])(0);
    }
    oldpts[num - 1] = (poly[0])(0);

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
    int count = 0;

    Real dist;
    for (int i = 0; i < num - 1; i++)
    {
        dist = norm(pts[count + 1] - pts[count]);
        if (dist <= (chdLenRange.hi())[0])
        {
            ++it;
            ++oit;
            count++;
            continue;
        }
        Vector<Real> chdt;
        split(knots[i], knots[i + 1], chdt);
        for (int j = 0; j < (int)chdt.size(); j++)
        {
            count++;
            ids.push_back(count);
        }
        count++;
        for (int j = chdt.size() - 1; j >= 0; j--)
        {
            Point opt = crv(chdt[j]);
            oit = oldpts.emplace(oit, opt);
            it = pts.emplace(it, TI->timeStep(v, opt, tn, dt));
        }
        it += chdt.size() + 1;
        oit += chdt.size() + 1;
    }
    //crv = fitCurve<Order>(oldpts, true);
    return ids;
}

template <int Dim, int Order>
void MARS<Dim, Order>::timeStep(const IT_VectorFunction<Dim> &v, YS &ys, Real tn, Real dt)
{
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    int id = 1;
    for (auto &crv : vcrv)
    {
        assert(crv.isClosed(tol));

        //pts: next time's points; oldpts: now time's points
        Vector<Point> oldpts;
        Vector<Polynomial<Order, Point>> poly = crv.getPolys();
        oldpts.resize(poly.size() + 1);
        for (int i = 0; i < (int)poly.size(); i++)
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

        int num = pts.size();
        Vector<Real> dist(num - 1);
        //dist.resize(num - 1);
        for (int i = 0; i < num - 1; i++)
        {
            dist[i] = norm(pts[i + 1] - pts[i]);
        }
        auto maxp = std::max_element(dist.begin(), dist.end());
        auto minp = std::min_element(dist.begin(), dist.end());

        std::cout << "Curve " << id << ":  Add " << splitids.size() << " Points, remove " << removeids.size() << " Points. "
                  << "Max chdlength: " << *maxp << "   Min chdlength: " << *minp << std::endl;
        id++;
    }
    ys = YS(SegmentedRealizableSpadjor<Order>(vcrv), tol);
    return;
}

template <int Dim, int Order>
void MARS<Dim, Order>::trackInterface(const IT_VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime, bool output, std::string fName, int opstride)
{
    std::cout << "chdLenRange: [" << (chdLenRange.lo())[0] << ", " << (chdLenRange.hi())[0] << "]" << std::endl
              << std::endl;
    Real T = StartTime;
    Real t = dt;
    int step = 1;
    while (T < EndTime)
    {
        if (EndTime - T < dt)
        {
            t = EndTime - T;
        }
        std::cout << "Step: " << step << "     timestep: " << t << std::endl;
        timeStep(v, ys, T, t);
        std::cout << std::endl;
        T += t;

        step++;
    }
    return;
}

template class MARS<2, 4>;
*/

#endif