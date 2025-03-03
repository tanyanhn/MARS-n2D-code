#ifndef _VECTORFORSURFACEDIFFUSIONFLOW_H_
#define _VECTORFORSURFACEDIFFUSIONFLOW_H_

#include "FlowHelper.h"
#include "VectorOnHypersurface.h"

template <int Dim, int Order>
class VectorForSurfaceDiffusionFlow : public VectorOnHypersurface<Dim> {
 private:
  template <class T>
  using Vector = std::vector<T>;

  using Point = Vec<Real, Dim>;

 public:
  const Vector<Point> operator()(const Vector<Point>& pts, Real t = 0) const;

  const Tensor<Real, 2> getJacobi(const Vector<Point>& pts, Real t = 0) const;
};

template <int Order>
class VectorForSurfaceDiffusionFlow<2, Order> : public VectorOnHypersurface<2> {
 private:
  template <class T>
  using Vector = std::vector<T>;

  using Point = Vec<Real, 2>;

 public:
  const Vector<Point> operator()(const Vector<Point>& pts, Real t = 0) const;

  const Tensor<Real, 2> getJacobi(const Vector<Point>& pts, Real t = 0) const;
};

template <int Order>
const Vector<Point> VectorForSurfaceDiffusionFlow<2, Order>::operator()(
    const Vector<Point>& pts, Real t) const {
  const int num = pts.size();
  const Curve<2, 4> crv = fitCurve<4>(pts, Curve<2, 4>::periodic);
  Vector<Point> der1 = calDer<Order + 2>(pts, crv, 1);
  Vector<Point> der2 = calDer<Order + 2>(pts, crv, 2);
  Vector<Real> kappa = Vector<Real>(num);
  for (int i = 0; i < num; i++)
    kappa[i] = -der2[i][0] * der1[i][1] + der2[i][1] * der1[i][0];
  Vector<Real> d2kappa = calDer<Order + 2>(kappa, crv, 2);
  Vector<Point> res = Vector<Point>(num);
  for (int i = 0; i < num; i++) {
    res[i][0] = d2kappa[i] * der1[i][1];
    res[i][1] = -d2kappa[i] * der1[i][0];
  }
  return res;
}

template <int Order>
const Tensor<Real, 2> VectorForSurfaceDiffusionFlow<2, Order>::getJacobi(
    const Vector<Point>& pts, Real t) const {
  const int num = pts.size();
  const Curve<2, 4> crv = fitCurve<4>(pts, Curve<2, 4>::periodic);
  Vector<Point> der1 = calDer<Order + 2>(pts, crv, 1);
  Vector<Point> der2 = calDer<Order + 2>(pts, crv, 2);
  Vector<Real> kappa = Vector<Real>(num);
  for (int i = 0; i < num; i++)
    kappa[i] = -der2[i][0] * der1[i][1] + der2[i][1] * der1[i][0];
  Vector<Real> d2kappa = calDer<Order + 2>(kappa, crv, 2);
  Vector<Real> arclength = calArcLength<Order / 2 + 2>(crv);
  Vector<Tensor<Real, 2> > dkappadx(num - 1);
  for (int i = 0; i < num - 1; i++)
    dkappadx[i] = calLocaldkappadx<Order + 2>(pts, crv, i, der1, der2);
  Tensor<Real, 2> res(Vec<int, 2>{2 * (num - 1), 2 * (num - 1)});
  for (int i = 0; i < 2 * num - 2; i++)
    for (int j = 0; j < 2 * num - 2; j++) res(i, j) = 0.0;
  for (int i = 0; i < num - 1; i++) {
    Vector<Real> Localcoes2 = calLocalFDcoes<Order + 2>(arclength, i, 2);
    Tensor<Real, 2> Localdderdx1 = calLocaldderdx1<Order + 2>(pts, crv, i);
    for (int m = i - Order / 2 - 1; m <= i + Order / 2 + 2; m++) {
      int tmpm = 0;
      if (m >= 0 && m < num - 1)
        tmpm = m;
      else if (m < 0)
        tmpm = num - 1 + m;
      else
        tmpm = m - num + 1;
      double tmp1 = Localcoes2[m - i + Order / 2 + 1] * der1[i][1];
      double tmp2 = Localcoes2[m - i + Order / 2 + 1] * (-der1[i][0]);
      for (int j = m - Order / 2 - 1; j <= m + Order / 2 + 2; j++) {
        int tmpj = 0;
        if (j >= 0 && j < num - 1)
          tmpj = j;
        else if (j < 0)
          tmpj = num - 1 + j;
        else
          tmpj = j - num + 1;
        res(i, tmpj) += tmp1 * dkappadx[tmpm](0, j - m + Order / 2 + 1);
        res(i, tmpj + num - 1) +=
            tmp1 * dkappadx[tmpm](1, j - m + Order / 2 + 1);
        res(i + num - 1, tmpj) +=
            tmp2 * dkappadx[tmpm](0, j - m + Order / 2 + 1);
        res(i + num - 1, tmpj + num - 1) +=
            tmp2 * dkappadx[tmpm](1, j - m + Order / 2 + 1);
      }
      const double Locald2kappa = d2kappa[i];
      res(i, tmpm) += Locald2kappa * Localdderdx1(2, m - i + Order / 2 + 1);
      res(i, tmpm + num - 1) +=
          Locald2kappa * Localdderdx1(3, m - i + Order / 2 + 1);
      res(i + num - 1, tmpm) +=
          Locald2kappa * (-Localdderdx1(0, m - i + Order / 2 + 1));
      res(i + num - 1, tmpm + num - 1) +=
          Locald2kappa * (-Localdderdx1(1, m - i + Order / 2 + 1));
    }
  }
  return res;
}

#endif
