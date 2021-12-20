#ifndef _EXPLICIT_RK_H_
#define _EXPLICIT_RK_H_

#include <vector>
#include "TimeIntegrator.h"
#include "RKButcher.h"

/// The compile time constants of numbers of stages

template <int Dim, RK::Type_Minor Type>
class ERK : public TimeIntegrator<Dim>
{

    template <class T>
    using Vector = std::vector<T>;

    using Point = Vec<Real, Dim>;

    using ButcherTab = ButcherTableau<RK::ERK, Type>;

    //enum{order = ButcherTab::nStages};

public:
    const int order = ButcherTab::order;

    void timeStep(const VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real dt)
    {
        int num = pts.size();
        Vector<Vector<Point>> step;
        step.resize(ButcherTab::nStages);
        Vector<Point> tmpts;
        Vector<Point> result = pts;

        for (int i = 0; i < ButcherTab::nStages; i++)
        {
            tmpts = pts;
            for (int j = 0; j < i; j++)
            {
                for (int k = 0; k < num; k++)
                {
                    tmpts[k] = tmpts[k] + step[j][k] * ButcherTab::a[i][j] * dt;
                }
            }
            step[i] = v(tmpts, tn + ButcherTab::c[i] * dt);
            for (int k = 0; k < num; k++)
            {
                result[k] = result[k] + step[i][k] * ButcherTab::b[i] * dt;
            }
        }
        pts = result;
    }
};

#endif