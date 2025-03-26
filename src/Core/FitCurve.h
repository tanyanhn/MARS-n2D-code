#pragma once
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "Core/Curve.h"

Eigen::VectorXd solveSystem(int n, const Eigen::SparseMatrix<double>& A,
                            const Eigen::VectorXd& rhs, Real tol = newtonTol(),
                            int maxDirectSize = 2000, int maxIter = 1000);

#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
inline Curve<2, 4>
fitCurveEigen(const std::vector<Vec<Real, 2>>& points,
              typename Curve<2, 4>::BCType type, Real tol = newtonTol()) {
  using Eigen::VectorXd;
  using std::array;
  using std::vector;
  using Triplet = Eigen::Triplet<double>;

  const int n = points.size();
  if (n < 2) return {};

  // Calculate cumulative distances as knots
  vector<Real> knots(n);
  knots[0] = 0;
  for (int i = 1; i < n; ++i) {
    Real dx = points[i][0] - points[i - 1][0];
    Real dy = points[i][1] - points[i - 1][1];
    knots[i] = knots[i - 1] + hypot(dx, dy);
  }

  // vector<array<array<Real, 4>, 2>> result(n - 1);
  std::vector<Polynomial<4, Vec<Real, 2>>> result(n - 1);
  // const Real total_length = knots.back();

  // Build and solve spline for each coordinate
  for (int coord = 0; coord < 2; ++coord) {
    VectorXd y(n);
    for (int i = 0; i < n; ++i) y[i] = points[i][coord];

    VectorXd h(n - 1);
    for (int i = 0; i < n - 1; ++i) h[i] = knots[i + 1] - knots[i];

    if (type == Curve<2, 4>::notAKnot) {
      if (n <= 3)
        throw std::runtime_error("Not-a-Knot condition not supported");
      // Not-a-Knot condition
      VectorXd rhs(n);
      Eigen::SparseMatrix<double> A(n, n);
      std::vector<Triplet> triplets;
      VectorXd alpha(n);
      VectorXd beta(n);
      VectorXd gamma(n);

      // 内部节点方程 (i=1 到 n-2)
      for (int i = 1; i < n - 1; ++i) {
        double hi_prev = h[i - 1];
        double hi = h[i];
        triplets.emplace_back(i, i - 1, hi_prev);
        triplets.emplace_back(i, i, 2 * (hi_prev + hi));
        triplets.emplace_back(i, i + 1, hi);
        rhs[i] =
            6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
      }

      // 第一个 Not-a-Knot 条件 (i=0)
      triplets.emplace_back(0, 0, h[1]);
      triplets.emplace_back(0, 1, -(h[0] + h[1]));
      triplets.emplace_back(0, 2, h[0]);
      rhs[0] = 0.0;

      // 最后一个 Not-a-Knot 条件 (i=n-1)
      triplets.emplace_back(n - 1, n - 3, h[n - 2]);
      triplets.emplace_back(n - 1, n - 2, -(h[n - 3] + h[n - 2]));
      triplets.emplace_back(n - 1, n - 1, h[n - 3]);
      rhs[n - 1] = 0.0;

      A.setFromTriplets(triplets.begin(), triplets.end());

      VectorXd M = solveSystem(n, A, rhs, tol);

      // Calculate coefficients
      for (int i = 0; i < n - 1; ++i) {
        const Real t1 = M[i + 1] - M[i];
        result[i][0][coord] = y[i];
        result[i][1][coord] =
            (y[i + 1] - y[i]) / h[i] - h[i] * (2 * M[i] + M[i + 1]) / 6;
        result[i][2][coord] = M[i] / 2;
        result[i][3][coord] = t1 / (6 * h[i]);
      }
    } else if (type == Curve<2, 4>::periodic) {
      if (VecCompare<Real, 2>(distTol()).compare(points[0], points[n - 1]) != 0)
        throw std::runtime_error("Periodic condition not satisfied");
      // Periodic condition
      int m = n - 1;
      VectorXd rhs(m);
      Eigen::SparseMatrix<double> A(m, m);
      std::vector<Triplet> triplets;

      // 内部节点方程 (i=1 到 n-2)
      for (int i = 1; i < m - 1; ++i) {
        double hi_prev = h[i - 1];
        double hi = h[i];

        // 二阶导数连续条件
        triplets.emplace_back(i, i - 1, hi_prev);
        triplets.emplace_back(i, i, 2 * (hi_prev + hi));
        triplets.emplace_back(i, i + 1, hi);
        rhs[i] = 6 * ((y[i + 1] - y[i]) / hi - (y[i] - y[i - 1]) / hi_prev);
      }

      // 周期性边界方程
      // 条件1: 二阶导数周期连续性（连接 S_0 和 S_{m-1}）
      triplets.emplace_back(0, m - 1, h[m - 1]);           // M_{m-1} 的系数
      triplets.emplace_back(0, 0, 2 * (h[m - 1] + h[0]));  // M_{0} 的系数
      triplets.emplace_back(0, 1, h[0]);                   // M_{1} 的系数
      rhs[0] = 6 * ((y[1] - y[0]) / h[0] - (y[0] - y[m - 1]) / h[m - 1]);

      triplets.emplace_back(m - 1, m - 2, h[m - 2]);  // M_{m-2} 的系数
      triplets.emplace_back(m - 1, m - 1,
                            2 * (h[m - 2] + h[m - 1]));  // M_{m-1} 的系数
      triplets.emplace_back(m - 1, 0, h[m - 1]);         // M_{0} 的系数
      rhs[m - 1] =
          6 * ((y[0] - y[m - 1]) / h[m - 1] - (y[m - 1] - y[m - 2]) / h[m - 2]);

      A.setFromTriplets(triplets.begin(), triplets.end());

      // VectorXd M = solveTridiagonal(alpha, beta, gamma, rhs);
      auto M = solveSystem(m, A, rhs, tol);

      // Calculate coefficients
      for (int i = 0; i < m; ++i) {
        const int next = (i + 1) % m;
        const Real t1 = M[next] - M[i];
        result[i][0][coord] = y[i];
        result[i][1][coord] =
            (y[next] - y[i]) / h[i] - h[i] * (2 * M[i] + M[next]) / 6;
        result[i][2][coord] = M[i] / 2;
        result[i][3][coord] = t1 / (6 * h[i]);
      }
    }
  }

  return Curve<2, 4>(std::move(knots), std::move(result));
}

inline Eigen::VectorXd solveSystem(int n, const Eigen::SparseMatrix<double>& A,
                                   const Eigen::VectorXd& rhs, Real tol,
                                   int maxDirectSize, int maxIter) {
  using namespace Eigen;

  // 智能阈值判定：基于矩阵非零元素数量
  const size_t nonZeros = A.nonZeros();
  const bool useDirect =
      (n <= maxDirectSize) && (nonZeros < 1e7);  // 防止超大稀疏矩阵用直接法

  if (useDirect) {
    // 小规模问题：使用带行列优化的稀疏LU
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;

    solver.compute(A);
    if (solver.info() == Success) {
      VectorXd result = solver.solve(rhs);
      if (solver.info() == Success) {
        return result;
      }
    }
    // 直接法失败时转用迭代法
  }

  // 中大规模问题：使用带预处理的迭代法
  BiCGSTAB<SparseMatrix<double>, DiagonalPreconditioner<double>> solver;

  // 自适应参数配置
  solver.setTolerance(tol);
  solver.setMaxIterations(maxIter);

  solver.compute(A);
  if (solver.info() != Success) {
    throw std::runtime_error("Iterative solver setup failed");
  }

  VectorXd result = solver.solve(rhs);
  // auto err = rhs - A * result;
  // std::cout << "err:" << err.norm() << std::endl;

  if (solver.info() != Success) {
    throw std::runtime_error("All iterative methods failed");
  }

  return result;
}
