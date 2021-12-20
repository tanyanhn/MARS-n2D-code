#ifndef _DIRK_H_
#define _DIRK_H_

#include <vector>
#include "TimeIntegrator.h"
#include "RKButcher.h"
#include "Core/Tensor.h"
#include "Core/TensorExpr.h"
#include "Core/TensorSlice.h"
#include "Core/Wrapper_LAPACKE.h"

template <int Dim>
std::vector<Vec<Real, Dim>> solveNewton(const VectorFunction<Dim> &v, Real coef, const Tensor<Real, 2> &Jacobi, const std::vector<Vec<Real, Dim>> &rhs, const std::vector<Vec<Real, Dim>> &pts, Real tol)
{
    using Point = Vec<Real, Dim>;

    int num = pts.size() - 1;
    auto res = pts;
    assert(2 * num == Jacobi.size()[0]);

    Tensor<Real, 2> mat = -Jacobi * coef;
    for (int i = 0; i < 2 * num; i++)
    {
        mat(i, i) += 1;
    }

    std::vector<Point> vel;
    Tensor<Real, 1> err(2 * num);
    Tensor<int, 1> ipiv(2 * num);

    vel = v(res);
    for (int j = 0; j < num; j++)
    {
        err(j) = res[j][0] - rhs[j][0] - coef * vel[j][0];
        err(j + num) = res[j][1] - rhs[j][1] - coef * vel[j][1];
    }

    for (int i = 0; i < 10; i++)
    {
        if (norm(err) < tol)
            break;
        auto info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 2 * num, 1, mat.data(), 2 * num, ipiv.data(), err.data(), 2 * num);
        for (int j = 0; j < num; j++)
        {
            res[j][0] -= err(j);
            res[j][1] -= err(j + num);
        }
        res[num][0] -= err(0);
        res[num][1] -= err(num);
        vel = v(res);
        for (int j = 0; j < num; j++)
        {
            err(j) = res[j][0] - rhs[j][0] - coef * vel[j][0];
            err(j + num) = res[j][1] - rhs[j][1] - coef * vel[j][1];
        }
    }
    return res;
}

template <int Dim, RK::Type_Minor Type>
class DIRK : public TimeIntegrator<Dim>
{

    template <class T>
    using Vector = std::vector<T>;

    using Point = Vec<Real, Dim>;

    using ButcherTab = ButcherTableau<RK::DIRK, Type>;

public:
    const int order = ButcherTab::order;

    void timeStep(const VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real dt)
    {
        int num = pts.size();
        Tensor<Real, 2> Jacobi = v.getJacobi(pts, tn);
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
            if (ButcherTab::a[i][i] != 0)
            {
                auto tmpts = solveNewton(v, ButcherTab::a[i][i] * dt, Jacobi, tmpts, pts, 1e-10);
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