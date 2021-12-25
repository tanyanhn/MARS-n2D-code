#ifndef _DIRK_H_
#define _DIRK_H_

#include <vector>
#include "TimeIntegrator.h"
#include "RKButcher.h"
#include "Core/Tensor.h"
#include "Core/TensorExpr.h"
#include "Core/TensorSlice.h"
#include "Core/Wrapper_LAPACKE.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

template <int Dim>
std::vector<Vec<Real, Dim>> solveNewton(const VectorFunction<Dim> &v, Real tn, Real coef, const Eigen::SparseMatrix<Real> &Jacobi, const std::vector<Vec<Real, Dim>> &rhs, const std::vector<Vec<Real, Dim>> &pts, Real tol)
{
    using Point = Vec<Real, Dim>;

    using SpMat = Eigen::SparseMatrix<Real>;

    int num = pts.size() - 1;
    auto res = pts;

    SpMat mat = -coef * Jacobi;
    for (int i = 0; i < Dim * num; i++)
    {
        mat.coeffRef(i, i) += 1;
    }

    Eigen::SparseLU<SpMat> lu(mat);
    
    //std::cout << mat << std::endl;

    std::vector<Point> vel;
    Eigen::VectorXd err(Dim * num);
    //Tensor<int, 1> ipiv(Dim * num);

    auto normXd = [](const Eigen::VectorXd err, int n) -> Real
    {
        Real res = 0;
        for(int i=0; i<n; i++)
        {
            res += pow(err(i),2);
        }
        return sqrt(res);
    };

    vel = v(res, tn);
    for (int j = 0; j < num; j++)
    {
        err(j) = res[j][0] - rhs[j][0] - coef * vel[j][0];
        err(j + num) = res[j][1] - rhs[j][1] - coef * vel[j][1];
        if (Dim == 3)
            err(j + 2 * num) = res[j][2] - rhs[j][2] - coef * vel[j][2];
    }

    for (int i = 0; i < 10; i++)
    {
        if (normXd(err, Dim * num) < tol)
            break;
        //auto info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 2 * num, 1, mat.data(), 2 * num, ipiv.data(), err.data(), 2 * num);
        auto x = lu.solve(err);
        if (lu.info() != Eigen::Success)
        {
            //std::cout << info << std::endl;
            throw std::runtime_error("solveNewton() - EigenSparseLU");
        }
        for (int j = 0; j < num; j++)
        {
            res[j][0] -= x(j);
            res[j][1] -= x(j + num);
            if (Dim == 3)
            {
                res[j][2] -= x(j + 2 * num);
            }
        }
        res[num][0] -= x(0);
        res[num][1] -= x(num);
        if (Dim == 3)
        {
            res[num][2] -= x(2 * num);
        }
        vel = v(res, tn);
        for (int j = 0; j < num; j++)
        {
            err(j) = res[j][0] - rhs[j][0] - coef * vel[j][0];
            err(j + num) = res[j][1] - rhs[j][1] - coef * vel[j][1];
            if (Dim == 3)
                err(j + 2 * num) = res[j][2] - rhs[j][2] - coef * vel[j][2];
        }
        if (i == 10)
        {
            exit(1);
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

    using SpMat = Eigen::SparseMatrix<Real>;

public:
    const int order = ButcherTab::order;

    void timeStep(const VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real dt)
    {
        int num = pts.size();
        //Tensor<Real, 2> jacobi = v.getJacobi(pts, tn);
        SpMat jacobi = v.getJacobi(pts, tn);
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
                tmpts = solveNewton<Dim>(v, tn + ButcherTab::c[i] * dt, ButcherTab::a[i][i] * dt, jacobi, tmpts, tmpts, 1e-10);
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