#ifndef _EXPLICIT_RK_H_
#define _EXPLICIT_RK_H_

#include <vector>
#include "TimeIntegrator.h"
#include "RKButcher.h"

/// The compile time constants of numbers of stages

template <int Dim, RK_Category2 Type>
class ExplicitRungeKutta : public TimeIntegrator<Dim>
{

    template <class T>
    using Vector = std::vector<T>;

    using Point = Vec<Real, Dim>;

    using ButcherTab = ButcherTableau<ERK, Type>;

public:

    virtual const Point timeStep(const VectorFunction<Dim> &v, const Point &pt, Real tn, Real dt)
    {
        Vector<Point> step;
        step.resize(ButcherTab::nStages);
        Point tmp;
        Point result(pt);
        for (int i = 0; i < ButcherTab::nStages; i++)
        {
            tmp = pt;
            for (int j = 0; j < i; j++)
            {
                tmp = tmp + step[j] * ButcherTab::a[i][j] * dt;
            }
            step[i] = v(tmp, tn + ButcherTab::c[i] * dt);
            result = result + step[i] * ButcherTab::b[i] * dt;
        }
        return result;
    }
};

#endif