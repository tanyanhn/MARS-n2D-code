#ifndef _IT_RUNGEKUTTABASE_
#define _IT_RUNGEKUTTABASE_

#include <vector>
#include "Core/Config.h"
#include "Core/Vec.h"
#include "TimeIntegrator.h"
#include "RKButcher.h"

/// The compile time constants of numbers of stages

template <int Dim, class RKButcher>
class ExplicitRungeKutta : public IT_TimeIntegrator<Dim>
{

    using RKC = RKButcher;

    template <class T>
    using Vector = std::vector<T>;

    using Point = Vec<Real, Dim>;

private:
    //static constexpr int nStages = RKC::nS;
    //Vector<Vector<Real>> a;
    //Vector<Real> b, c;

public:
    /*ExplicitRungeKutta()
    {
        a.resize(nStages);
        for (int i = 0; i < nStages; i++)
        {
            a[i].resize(nStages);
        }
        b.resize(nStages);
        c.resize(nStages);
        for (int i = 0; i < nStages; i++)
        {
            for (int j = 0; j < nStages; j++)
                a[i][j] = static_cast<Real>(RKC::a[i][j]) / RKC::a[i][nStages];
            c[i] = static_cast<Real>(RKC::c[i]) / RKC::c[nStages];
        }
        for (int j = 0; j < nStages; j++)
        {
            b[j] = static_cast<Real>(RKC::a[nStages][j]) / RKC::a[nStages][nStages];
        }
    }*/

    virtual const Point timeStep(const IT_VectorFunction<Dim> &v, const Point &pt, Real tn, Real dt)
    {
        Vector<Point> step;
        step.resize(RKC::nS);
        Point tmp;
        Point result(pt);
        for (int i = 0; i < RKC::nS; i++)
        {
            tmp = pt;
            for (int j = 0; j < i; j++)
            {
                tmp = tmp + step[j] * RKC::a[i][j] / RKC::a[i][RKC::nS] * dt;
            }
            step[i] = v(tmp, tn + RKC::c[i] / RKC::c[RKC::nS] * dt);
            result = result + step[i] * RKC::b[i] / RKC::b[RKC::nS] * dt;
        }
        return result;
    }
};

#endif