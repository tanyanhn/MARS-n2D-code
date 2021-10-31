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

Real tol = 1e-15;

template <int Dim, int order>
class MARS
{
    using YS = YinSet<Dim, order>;

    template <class T>
    using Vector = std::vector<T>;

    using rVec = Vec<Real,Dim>;

    using Crv = Curve<Dim, order>;

private:
    IT_TimeIntegrator<Dim> *TI;
    Interval<1> chdLenRange;

public:
    MARS() = delete;
    MARS(IT_TimeIntegrator<Dim> *_TI, Real hL, Real rtiny = 0.1):TI(_TI),chdLenRange(Vec<Real,1>(rtiny*hL),Vec<Real,1>(hL))
    {}

    bool removeSmallEdges(Vector<rVec> &pts, Vector<rVec> &oldpts);

    bool splitLongEdges(const IT_VectorFunction<Dim> &v, Vector<rVec> &pts, Vector<rVec> &oldpts, const Curve<Dim, order> &crv, Real TN, Real dt);

    void discreteFlowMap(const IT_VectorFunction<Dim> &v, YS &ys, Real TN, Real dt);

    void trackInterface(const IT_VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime);
};

template <int Dim, int order>
bool MARS<Dim, order>::removeSmallEdges(Vector<rVec> &pts, Vector<rVec> &oldpts)
{
    assert(pts.size() == oldpts.size());
    size_t Num = pts.size();
    Vector<Real> dist;
    dist.resize(Num - 1);
    for (size_t i = 0; i < Num - 1; i++)
    {
        dist[i] = norm(pts[i + 1] - pts[i]);
    }

    //mark weather the pre-node is candidate to be erased.
    bool predelete = false; 
    auto it = pts.begin();
    auto oit = oldpts.begin();
    ++oit;
    ++it;
    size_t i = 1;
    if (dist[0] < (chdLenRange.lo())[0])
    {
        predelete = true;
        it = pts.erase(it);
        oit = oldpts.erase(oit);
        i++;
    }

    while (i < Num - 2)
    {
        if (dist[i] >= (chdLenRange.lo())[0])
        {
            predelete = false;
            ++it;
            ++oit;
            i++;
        }
        else
        {
            if (dist[i - 1] <= dist[i + 1])
            {
                if (predelete == false)
                {
                    it = pts.erase(it);
                    oit = oldpts.erase(oit);
                    i++;
                    predelete = true;
                }
                else
                {
                    ++it;
                    ++oit;
                    i++;
                    predelete = false;
                }
            }
            else
            {
                it = pts.erase(it + 1);
                oit = oldpts.erase(oit + 1);
                i += 2;
                predelete = true;
            }
        }
    }
    if (i == Num - 1)
    {
        return Num != pts.size();
    }
    else
    {
        if (predelete == true || dist[i] >= (chdLenRange.lo())[0])
            return Num != pts.size();
        it = pts.erase(it);
        oit = oldpts.erase(oit);
        return true;
    }
}

template <int Dim, int order>
bool MARS<Dim, order>::splitLongEdges(const IT_VectorFunction<Dim> &v, Vector<rVec> &pts, Vector<rVec> &oldpts, const Curve<Dim, order> &crv, Real TN, Real dt)
{
    assert(pts.size() == oldpts.size());
    Vector<Polynomial<4, rVec>> poly = crv.getPolys();
    Vector<Real> knots = crv.getKnots();
    size_t Num = pts.size();
    Vector<Real> dist;
    dist.resize(Num - 1);
    auto it = pts.begin();
    auto oit = oldpts.begin();
    for (size_t i = 0; i < Num - 1; i++)
    {
        dist[i] = norm(pts[i + 1] - pts[i]);
    }
    for (size_t i = 0; i < Num - 1; i++)
    {
        ++it;
        ++oit;
        if (dist[i] > (chdLenRange.hi())[0])
        {
            int N = (int)ceil(dist[i] / (chdLenRange.hi())[0]);
            Real dl = (knots[i + 1] - knots[i]) / N;
            for (int j = 1; j < N; j++)
            {
                rVec prept = (poly[i])(dl * j);
                oit = oldpts.emplace(oit, prept);
                ++oit;
                TI->timeStep(v, prept, TN, dt);
                it = pts.emplace(it, prept);
                ++it;
            }
        }
    }
    return Num != pts.size();
}

template <int Dim, int order>
void MARS<Dim, order>::discreteFlowMap(const IT_VectorFunction<Dim> &v, YS &ys, Real NowTime, Real dt)
{
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    for (auto &crv : vcrv)
    {
        assert(crv.isClosed(tol));

        //pts: next time's points; oldpts: now time's points
        Vector<rVec> oldpts;
        Vector<Polynomial<4, rVec>> poly = crv.getPolys();
        oldpts.resize(poly.size() + 1);
        for (size_t i = 0; i < poly.size(); i++)
        {
            oldpts[i] = (poly[i])[0];
        }
        oldpts[poly.size()] = (poly[0])[0];
        Vector<rVec> pts = oldpts;
        TI->multyTimeStep(v, pts, NowTime, dt);

        //get the next time's points
        

        //remove small edges, update pts and old pts simultaneously
        for (int i = 0; i < 3; i++)
        {
            if (!removeSmallEdges(pts, oldpts))
                break;
        }
        
        //regenerate the now time's curve
        crv = fitCurve<4>(oldpts, true);

        //split long edges, 
        for (int i = 0; i < 3; i++)
        {
            if (!splitLongEdges(v, pts, oldpts, crv, NowTime, dt))
                break;
            crv = fitCurve<4>(oldpts, true);
        }

        crv = fitCurve<4>(pts, true);
    }
    ys = YS(SegmentedRealizableSpadjor<4>(vcrv),tol);
    return;
}

template <int Dim, int order>
void MARS<Dim, order>::trackInterface(const IT_VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime)
{
    Real T = StartTime;
    Real t = dt;
    while (T < EndTime)
    {
        if(EndTime-T<dt)
        {
            t = EndTime - T;
        }
        discreteFlowMap(v, ys, T, t);
        T += t;
    }
    return;
}

template class MARS<2,4>;

#endif