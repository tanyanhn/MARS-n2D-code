#include "FlowHelper.h"

#include <cmath>

#include "Core/numlib.h"

template <int Order>
Real calArcLength(const Polynomial<4, Point>& poly, Real start, Real end) {
  auto dxdt = getComp(poly, 0).der();
  auto dydt = getComp(poly, 1).der();
  auto ndsdt = [&](Real t) {
    Vec<Real, 2> ds{dxdt(t), dydt(t)};
    return norm(ds, 2);
  };
  return quad<Order>(ndsdt, start, end);
}

template <int Order>
Vector<Real> calArcLength(const Crv& crv) {
  const Vector<Polynomial<4, Point> >& polys = crv.getPolys();
  const Vector<Real>& knots = crv.getKnots();
  const int num = knots.size();
  Vector<Real> res(num);
  res[0] = 0.0;
  for (int i = 1; i < num; i++)
    res[i] = res[i - 1] +
             calArcLength<Order>(polys[i - 1], 0, knots[i] - knots[i - 1]);
  return res;
}

Real adjustpoint(const Vector<Real>& arclength, const int j) {
  const int num = arclength.size();
  Real tmpj = 0.0;
  if (j >= 0 && j < num)
    tmpj = arclength[j];
  else if (j < 0)
    tmpj = arclength[j + num - 1] - arclength[num - 1];
  else
    tmpj = arclength[num - 1] + arclength[j - num + 1];
  return tmpj;
}

// template <int Order>
// Vector<Real> calLocalFD1coes(const Vector<Real>& arclength, const int
// i)
// {
//   assert(Order%2 == 0);
//   Vector<Real> res(Order+1);
//   for (int j = i - Order/2 ; j <= i + Order/2 ; j++){
//       Real tmpj = adjustpoint(arclength,j);
//       Polynomial<Order+1,Real> poly;
//       poly[0] = 1.0;
//       for (int k = i - Order/2 ; k <= i + Order/2 ; k++){
//         if (k != j){
//           Real tmpk = adjustpoint(arclength,k);
//           Polynomial<2,Real>
//             tmppoly{-(tmpk-arclength[i])/(tmpj-tmpk),1/(tmpj-tmpk)};
//           auto polyl = poly*tmppoly;
//           for (int m = 0 ; m < std::min(k - i + Order/2 + 2,Order+1); m++)
//             poly[m] = polyl[m];
//         }
//       }
//       auto dpoly = poly.der();
//       res[j-i+Order/2] = dpoly(0.0);
//     }
//   return res;
// }

// template <int Order>
// Vector<Real> calLocalFD2coes(const Vector<Real>& arclength, const int
// i)
// {
//   assert(Order%2 == 0);
//   Vector<Real> res(Order+2);
//   for (int j = i - Order/2 ; j <= i + Order/2 + 1 ; j++){
//     Real tmpj = adjustpoint(arclength,j);
//       Polynomial<Order+2,Real> poly;
//       poly[0] = 1.0;
//       for (int k = i - Order/2 ; k <= i + Order/2 + 1 ; k++){
//         if (k != j){
//           Real tmpk = adjustpoint(arclength,k);
//           Polynomial<2,Real>
//             tmppoly{-(tmpk-arclength[i])/(tmpj-tmpk),1/(tmpj-tmpk)};
//           auto polyl = poly*tmppoly;
//           for (int m = 0 ; m < std::min(k - i + Order/2 + 2,Order+2); m++)
//             poly[m] = polyl[m];
//         }
//       }
//       auto ddpoly = poly.der().der();
//       res[j-i+Order/2] = ddpoly(0.0);
//     }
//   return res;
// }

// template <int Order>
// Vector<Real> calLocalFD3coes(const Vector<Real>& arclength, const int
// i)
// {
//   assert(Order%2 == 0);
//   Vector<Real> res(Order+3);
//   for (int j = i - Order/2 - 1 ; j <= i + Order/2 + 1 ; j++){
//       Real tmpj = adjustpoint(arclength,j);
//       Polynomial<Order+3,Real> poly;
//       poly[0] = 1.0;
//       for (int k = i - Order/2 - 1 ; k <= i + Order/2 + 1 ; k++){
//         if (k != j){
//           Real tmpk = adjustpoint(arclength,k);
//           Polynomial<2,Real>
//             tmppoly{-(tmpk-arclength[i])/(tmpj-tmpk),1/(tmpj-tmpk)};
//           auto polyl = poly*tmppoly;
//           for (int m = 0 ; m < std::min(k - i + Order/2 + 3,Order+3); m++)
//             poly[m] = polyl[m];
//         }
//       }
//       auto ddpoly = poly.der().der().der();
//       res[j-i+Order/2+1] = ddpoly(0.0);
//     }
//   return res;
// }

// template <int Order>
// Vector<Real> calLocalFD4coes(const Vector<Real>& arclength, const int
// i)
// {
//   assert(Order%2 == 0);
//   Vector<Real> res(Order+4);
//   for (int j = i - Order/2 - 1 ; j <= i + Order/2 + 2 ; j++){
//       Real tmpj = adjustpoint(arclength,j);
//       Polynomial<Order+4,Real> poly;
//       poly[0] = 1.0;
//       for (int k = i - Order/2 -1 ; k <= i + Order/2 + 2 ; k++){
//         if (k != j){
//           Real tmpk = adjustpoint(arclength,k);
//           Polynomial<2,Real>
//             tmppoly{-(tmpk-arclength[i])/(tmpj-tmpk),1/(tmpj-tmpk)};
//           auto polyl = poly*tmppoly;
//           for (int m = 0 ; m < std::min(k - i + Order/2 + 4,Order+4); m++)
//             poly[m] = polyl[m];
//         }
//       }
//       auto ddpoly = poly.der().der().der().der();
//       res[j-i+Order/2+1] = ddpoly(0.0);
//     }
//   return res;
// }

template <int Order>
Vector<Real> calLocalFDcoes(const Vector<Real>& arclength, const int i,
                            const int M) {
  assert(Order % 2 == 0);
  Tensor<Real, 3> table(Vec<int, 3>{M + 1, Order + M, Order + M});
  table(0, 0, 0) = 1.0;
  double c1 = 1.0;
  for (int n = 1; n < Order + M; n++) {
    double c2 = 1.0;
    int j = i - Order / 2 - (M - 1) / 2 + n;
    Real tmpj = adjustpoint(arclength, j);
    for (int mu = 0; mu < n; mu++) {
      int k = i - Order / 2 - (M - 1) / 2 + mu;
      Real tmpk = adjustpoint(arclength, k);
      double c3 = tmpj - tmpk;
      c2 = c2 * c3;
      if (n <= M) table(n, n - 1, mu) = 0.0;
      for (int m = 0; m <= std::min(n, M); m++) {
        if (m > 0)
          table(m, n, mu) = ((tmpj - arclength[i]) * table(m, n - 1, mu) -
                             m * table(m - 1, n - 1, mu)) /
                            c3;
        else
          table(m, n, mu) = ((tmpj - arclength[i]) * table(m, n - 1, mu)) / c3;
      }
    }
    int jj = i - Order / 2 - (M - 1) / 2 + n - 1;
    Real tmpjj = adjustpoint(arclength, jj);
    for (int m = 0; m <= std::min(n, M); m++) {
      if (m > 0)
        table(m, n, n) = c1 / c2 *
                         (m * table(m - 1, n - 1, n - 1) -
                          (tmpjj - arclength[i]) * table(m, n - 1, n - 1));
      else
        table(m, n, n) =
            c1 / c2 * (-(tmpjj - arclength[i]) * table(m, n - 1, n - 1));
    }
    c1 = c2;
  }
  Vector<Real> res(Order + M);
  for (int j = 0; j < Order + M; j++) res[j] = table(M, Order + M - 1, j);
  return res;
}

template <int Order>
Tensor<Real, 2> calFDcoes(const Vector<Real>& arclength, const int k) {
  assert(Order % 2 == 0);
  const int num = arclength.size();
  Tensor<Real, 2> res(Vec<int, 2>{Order + k, num - 1});
  for (int i = 0; i < num - 1; i++) {
    Vector<Real> LocalFDcoes = calLocalFDcoes<Order>(arclength, i, k);
    for (int j = 0; j < Order + k; j++) res(j, i) = LocalFDcoes[j];
  }
  return res;
}

// template <int Order>
// Tensor<Real,2> calFD1coes(const Vector<Real>& arclength)
// {
//   assert(Order%2 == 0);
//   const int num = arclength.size();
//   Tensor<Real,2> res(Vec<int,2> {Order+1,num-1});
//   for (int i = 0 ; i < num - 1 ; i++){
//     Vector<Real> LocalFD1coes = calLocalFD1coes<Order>(arclength,i);
//     for (int j = 0 ; j < Order + 1 ; j++)
//       res(j,i)=LocalFD1coes[j];
//   }
//   return res;
// }

// template <int Order>
// Tensor<Real,2> calFD2coes(const Vector<Real>& arclength)
// {
//   assert(Order%2 == 0);
//   const int num = arclength.size();
//   Tensor<Real,2> res(Vec<int,2> {Order+2,num-1});
//   for (int i = 0 ; i < num - 1 ; i++){
//     Vector<Real> LocalFD2coes = calLocalFD2coes<Order>(arclength,i);
//     for (int j = 0 ; j < Order + 2 ; j++)
//       res(j,i)=LocalFD2coes[j];
//   }
//   return res;
// }

// template <int Order>
// Tensor<Real,2> calFD3coes(const Vector<Real>& arclength)
// {
//   assert(Order%2 == 0);
//   const int num = arclength.size();
//   Tensor<Real,2> res(Vec<int,2> {Order+3,num-1});
//   for (int i = 0 ; i < num - 1 ; i++){
//     Vector<Real> LocalFD3coes = calLocalFD3coes<Order>(arclength,i);
//     for (int j = 0 ; j < Order + 3 ; j++)
//       res(j,i)=LocalFD3coes[j];
//   }
//   return res;
// }

// template <int Order>
// Tensor<Real,2> calFD4coes(const Vector<Real>& arclength)
// {
//   assert(Order%2 == 0);
//   const int num = arclength.size();
//   Tensor<Real,2> res(Vec<int,2> {Order+4,num-1});
//   for (int i = 0 ; i < num - 1 ; i++){
//     Vector<Real> LocalFD4coes = calLocalFD4coes<Order>(arclength,i);
//     for (int j = 0 ; j < Order + 4 ; j++)
//       res(j,i)=LocalFD4coes[j];
//   }
//   return res;
// }

template <int Order, class T>
Vector<T> calDer(const Vector<T>& datas, const Crv& crv, const int k) {
  assert(Order % 2 == 0);
  Vector<Real> arclength = calArcLength<Order / 2 + 1>(crv);
  Tensor<Real, 2> FDcoes = calFDcoes<Order>(arclength, k);
  const int num = datas.size();
  Vector<T> res(num);
  const int numterms = (FDcoes.size())[0];
  for (int i = 0; i < num - 1; i++) {
    for (int j = 0; j < numterms; j++) {
      int tmp = 0;
      const int index = i - Order / 2 + j - (k - 1) / 2;
      if (index >= 0 && index < num)
        tmp = index;
      else if (index < 0)
        tmp = index + num - 1;
      else
        tmp = index - num + 1;
      res[i] = res[i] + datas[tmp] * FDcoes(j, i);
    }
  }
  res[num - 1] = res[0];
  return res;
}

Vector<Real> calLocaldsdx(const Point p1, const Point p2) {
  const double s = norm(p1 - p2, 2);
  Vector<Real> res(4);
  res[0] = (p1[0] - p2[0]) / s;
  res[1] = -res[0];
  res[2] = (p1[1] - p2[1]) / s;
  res[3] = -res[2];
  return res;
}

Tensor<Real, 2> caldsdx(const Vector<Point>& pts) {
  const int num = pts.size();
  Tensor<Real, 2> res(Vec<int, 2>{4, num - 1});
  for (int i = 0; i < num - 1; i++) {
    const double s = norm(pts[i] - pts[i + 1], 2);
    res(0, i) = (pts[i][0] - pts[i + 1][0]) / s;
    res(1, i) = -res(0, i);
    res(2, i) = (pts[i][1] - pts[i + 1][1]) / s;
    res(3, i) = -res(2, i);
  }
  return res;
}

Vector<Real> adjustarc(const Vector<Real>& arc, const int tmpl,
                       const double dx) {
  const int num = arc.size();
  Vector<Real> tmparc(num - 1);
  for (int i = 0; i < num - 1; i++) tmparc[i] = arc[i + 1] - arc[i];
  tmparc[tmpl] += dx;
  Vector<Real> res(num);
  res[0] = 0.0;
  for (int i = 1; i < num; i++) res[i] = res[i - 1] + tmparc[i - 1];
  return res;
}

template <int Order>
Tensor<Real, 2> calLocaldads1(const Vector<Point>& pts, const Crv& crv,
                              const int i, const double dx) {
  assert(Order % 2 == 0);
  Vector<Real> arclength = calArcLength<Order / 2 + 1>(crv);
  const int num = arclength.size();
  Tensor<Real, 2> res(Vec<int, 2>{Order + 1, Order});
  // Before small increment
  Vector<Real> coesBefore = calLocalFDcoes<Order>(arclength, i, 1);
  // After small increment
  for (int l = i - Order / 2; l <= i + Order / 2 - 1; l++) {
    int tmpl = 0;
    if (l >= 0 && l < num - 1)
      tmpl = l;
    else if (l < 0)
      tmpl = num - 1 + l;
    else
      tmpl = l - num + 1;
    Vector<Real> newarclength = adjustarc(arclength, tmpl, dx);
    Vector<Real> coesAfter = calLocalFDcoes<Order>(newarclength, i, 1);
    for (int j = i - Order / 2; j <= i + Order / 2; j++)
      res(j - i + Order / 2, l - i + Order / 2) =
          (coesAfter[j - i + Order / 2] - coesBefore[j - i + Order / 2]) / dx;
  }
  return res;
}

template <int Order>
Tensor<Real, 2> calLocaldads2(const Vector<Point>& pts, const Crv& crv,
                              const int i, const double dx) {
  assert(Order % 2 == 0);
  Vector<Real> arclength = calArcLength<Order / 2 + 1>(crv);
  const int num = arclength.size();
  Tensor<Real, 2> res(Vec<int, 2>{Order + 2, Order + 1});
  // Before small increment
  Vector<Real> coesBefore = calLocalFDcoes<Order>(arclength, i, 2);
  // After small increment
  for (int l = i - Order / 2; l <= i + Order / 2; l++) {
    int tmpl = 0;
    if (l >= 0 && l < num - 1)
      tmpl = l;
    else if (l < 0)
      tmpl = num - 1 + l;
    else
      tmpl = l - num + 1;
    Vector<Real> newarclength = adjustarc(arclength, tmpl, dx);
    Vector<Real> coesAfter = calLocalFDcoes<Order>(newarclength, i, 2);
    for (int j = i - Order / 2; j <= i + Order / 2 + 1; j++)
      res(j - i + Order / 2, l - i + Order / 2) =
          (coesAfter[j - i + Order / 2] - coesBefore[j - i + Order / 2]) / dx;
  }
  return res;
}

template <int Order>
Tensor<Real, 3> caldads1(const Vector<Point>& pts, const Crv& crv,
                         const double dx) {
  assert(Order % 2 == 0);
  const int num = pts.size();
  Tensor<Real, 3> res(Vec<int, 3>{num, Order + 1, Order});
  for (int i = 0; i < num; i++) {
    Tensor<Real, 2> tmptensor = calLocaldads1<Order>(pts, crv, i, dx);
    for (int j = 0; j < Order + 1; j++)
      for (int k = 0; k < Order; k++) res(i, j, k) = tmptensor(j, k);
  }
  return res;
}

template <int Order>
Tensor<Real, 3> caldads2(const Vector<Point>& pts, const Crv& crv,
                         const double dx) {
  assert(Order % 2 == 0);
  const int num = pts.size();
  Tensor<Real, 3> res(Vec<int, 3>{num, Order + 2, Order + 1});
  for (int i = 0; i < num; i++) {
    Tensor<Real, 2> tmptensor = calLocaldads2<Order>(pts, crv, i, dx);
    for (int j = 0; j < Order + 2; j++)
      for (int k = 0; k < Order + 1; k++) res(i, j, k) = tmptensor(j, k);
  }
  return res;
}

// template <int Order>
// Vector<Real> calLocaldderdx2(const Vector<Point>& pts, const Crv& crv,
//                              const int i, const int di){
//   assert(Order%2 == 0);
//   Vector<Real> res{0,0,0,0};
//   if (di < -Order/2 || di > Order/2+1 )
//     return res;
//   const int num = pts.size();
//   Vector<Real> arclength = calArcLength<Order/2+1>(crv);

//   // cal FD coes we need
//   const int j = i + di;
//   Real tmpj = 0;
//   if (j >= 0 && j < num)
//     tmpj = arclength[j];
//   else if (j < 0)
//     tmpj = arclength[j+num-1] - arclength[num-1];
//   else
//     tmpj = arclength[num-1] + arclength[j-num+1];
//   Polynomial<Order+2,Real> poly;
//   poly[0] = 1.0;
//   for (int k = i - Order/2 ; k <= i + Order/2 + 1 ; k++){
//     if (k != j){
//       Real tmpk = 0;
//       if (k >= 0 && k < num)
//         tmpk = arclength[k];
//       else if (k < 0)
//         tmpk = arclength[k+num-1] - arclength[num-1];
//       else
//         tmpk = arclength[num-1] + arclength[k-num+1];
//       Polynomial<2,Real>
//         tmppoly{-(tmpk-arclength[i])/(tmpj-tmpk),1/(tmpj-tmpk)};
//       auto polyl = poly*tmppoly;
//       for (int m = 0 ; m < std::min(k - i + Order/2 + 2,Order+2); m++)
//         poly[m] = polyl[m];
//     }
//   }
//   auto ddpoly = poly.der().der();
//   Real Localcoe = ddpoly(0.0);
//   // prepare for building results: dsdx
//   int indextmpj = 0;
//   if (j >= 0 && j < num - 1)
//     indextmpj = j;
//   else if (j < 0)
//     indextmpj = num-1+j;
//   else
//     indextmpj = j-num+1;
//   Vector<Real> Localdsdx1;
//   if (indextmpj != 0)
//     Localdsdx1 = calLocaldsdx(pts[indextmpj-1],pts[indextmpj]);
//   else
//     Localdsdx1 = calLocaldsdx(pts[num-2],pts[num-1]);
//   Vector<Real> Localdsdx2;
//   if (indextmpj != num-1)
//     Localdsdx2 = calLocaldsdx(pts[indextmpj],pts[indextmpj+1]);
//   else
//     Localdsdx2 = calLocaldsdx(pts[0],pts[1]);
//   // prepare for building results: dads
//   Tensor<Real,2> Localdads = calLocaldads2<Order>(pts,crv,i,1e-6);
//   // building results
//   res[0]+=Localcoe;
//   res[3]+=Localcoe;
//   for (int k = i - Order/2; k <= i + Order/2 + 1 ; k++){
//     int tmpk = 0;
//     if (k >= 0 && k < num - 1)
//       tmpk = k;
//     else if (k < 0)
//       tmpk = num-1+k;
//     else
//       tmpk = k-num+1;
//     if (j != i - Order/2){
//       res[0]+=pts[tmpk][0]*Localdads(k-i+Order/2,j-i+Order/2-1)*Localdsdx1[1];
//       res[1]+=pts[tmpk][0]*Localdads(k-i+Order/2,j-i+Order/2-1)*Localdsdx1[3];
//       res[2]+=pts[tmpk][1]*Localdads(k-i+Order/2,j-i+Order/2-1)*Localdsdx1[1];
//       res[3]+=pts[tmpk][1]*Localdads(k-i+Order/2,j-i+Order/2-1)*Localdsdx1[3];
//     }
//     if (j != i + Order/2 + 1){
//       res[0]+=pts[tmpk][0]*Localdads(k-i+Order/2,j-i+Order/2)*Localdsdx2[0];
//       res[1]+=pts[tmpk][0]*Localdads(k-i+Order/2,j-i+Order/2)*Localdsdx2[2];
//       res[2]+=pts[tmpk][1]*Localdads(k-i+Order/2,j-i+Order/2)*Localdsdx2[0];
//       res[3]+=pts[tmpk][1]*Localdads(k-i+Order/2,j-i+Order/2)*Localdsdx2[2];
//     }
//   }
//   return res;
// }

template <int Order>
Tensor<Real, 2> calLocaldderdx2(const Vector<Point>& pts, const Crv& crv,
                                const int i) {
  assert(Order % 2 == 0);
  Tensor<Real, 2> res(Vec<int, 2>{4, Order + 2});
  for (int l = 0; l < 4; l++)
    for (int m = 0; m < Order + 2; m++) res(l, m) = 0.0;
  const int num = pts.size();
  Vector<Real> arclength = calArcLength<Order / 2 + 1>(crv);
  // cal FD coes we need
  Vector<Real> coes = calLocalFDcoes<Order>(arclength, i, 2);
  // prepare for building results: dsdx
  Vector<Vector<Real> > Localdsdx1(Order + 2);
  Vector<Vector<Real> > Localdsdx2(Order + 2);
  for (int j = i - Order / 2; j <= i + Order / 2 + 1; j++) {
    int indextmpj = 0;
    if (j >= 0 && j < num - 1)
      indextmpj = j;
    else if (j < 0)
      indextmpj = num - 1 + j;
    else
      indextmpj = j - num + 1;
    if (indextmpj != 0)
      Localdsdx1[j - i + Order / 2] =
          calLocaldsdx(pts[indextmpj - 1], pts[indextmpj]);
    else
      Localdsdx1[j - i + Order / 2] = calLocaldsdx(pts[num - 2], pts[num - 1]);
    if (indextmpj != num - 1)
      Localdsdx2[j - i + Order / 2] =
          calLocaldsdx(pts[indextmpj], pts[indextmpj + 1]);
    else
      Localdsdx2[j - i + Order / 2] = calLocaldsdx(pts[0], pts[1]);
  }
  // prepare for building results: dads
  Tensor<Real, 2> Localdads = calLocaldads2<Order>(pts, crv, i, 1e-6);
  // building results
  for (int j = i - Order / 2; j <= i + Order / 2 + 1; j++) {
    res(0, j - i + Order / 2) += coes[j - i + Order / 2];
    res(3, j - i + Order / 2) += coes[j - i + Order / 2];
    for (int k = i - Order / 2; k <= i + Order / 2 + 1; k++) {
      int tmpk = 0;
      if (k >= 0 && k < num - 1)
        tmpk = k;
      else if (k < 0)
        tmpk = num - 1 + k;
      else
        tmpk = k - num + 1;
      if (j != i - Order / 2) {
        res(0, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][1];
        res(1, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][3];
        res(2, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][1];
        res(3, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][3];
      }
      if (j != i + Order / 2 + 1) {
        res(0, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][0];
        res(1, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][2];
        res(2, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][0];
        res(3, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][2];
      }
    }
  }
  return res;
}

template <int Order>
Tensor<Real, 2> calLocaldderdx1(const Vector<Point>& pts, const Crv& crv,
                                const int i) {
  assert(Order % 2 == 0);
  Tensor<Real, 2> res(Vec<int, 2>{4, Order + 1});
  for (int l = 0; l < 4; l++)
    for (int m = 0; m < Order + 1; m++) res(l, m) = 0.0;
  const int num = pts.size();
  Vector<Real> arclength = calArcLength<Order / 2 + 1>(crv);
  // cal FD coes we need
  Vector<Real> coes = calLocalFDcoes<Order>(arclength, i, 1);
  // prepare for building results: dsdx
  Vector<Vector<Real> > Localdsdx1(Order + 1);
  Vector<Vector<Real> > Localdsdx2(Order + 1);
  for (int j = i - Order / 2; j <= i + Order / 2; j++) {
    int indextmpj = 0;
    if (j >= 0 && j < num - 1)
      indextmpj = j;
    else if (j < 0)
      indextmpj = num - 1 + j;
    else
      indextmpj = j - num + 1;
    if (indextmpj != 0)
      Localdsdx1[j - i + Order / 2] =
          calLocaldsdx(pts[indextmpj - 1], pts[indextmpj]);
    else
      Localdsdx1[j - i + Order / 2] = calLocaldsdx(pts[num - 2], pts[num - 1]);
    if (indextmpj != num - 1)
      Localdsdx2[j - i + Order / 2] =
          calLocaldsdx(pts[indextmpj], pts[indextmpj + 1]);
    else
      Localdsdx2[j - i + Order / 2] = calLocaldsdx(pts[0], pts[1]);
  }
  // prepare for building results: dads
  Tensor<Real, 2> Localdads = calLocaldads1<Order>(pts, crv, i, 1e-6);
  // building results
  for (int j = i - Order / 2; j <= i + Order / 2; j++) {
    res(0, j - i + Order / 2) += coes[j - i + Order / 2];
    res(3, j - i + Order / 2) += coes[j - i + Order / 2];
    for (int k = i - Order / 2; k <= i + Order / 2; k++) {
      int tmpk = 0;
      if (k >= 0 && k < num - 1)
        tmpk = k;
      else if (k < 0)
        tmpk = num - 1 + k;
      else
        tmpk = k - num + 1;
      if (j != i - Order / 2) {
        res(0, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][1];
        res(1, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][3];
        res(2, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][1];
        res(3, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2 - 1) *
            Localdsdx1[j - i + Order / 2][3];
      }
      if (j != i + Order / 2) {
        res(0, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][0];
        res(1, j - i + Order / 2) +=
            pts[tmpk][0] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][2];
        res(2, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][0];
        res(3, j - i + Order / 2) +=
            pts[tmpk][1] * Localdads(k - i + Order / 2, j - i + Order / 2) *
            Localdsdx2[j - i + Order / 2][2];
      }
    }
  }
  return res;
}

template <int Order>
Tensor<Real, 2> calLocaldkappadx(const Vector<Point>& pts, const Crv& crv,
                                 const int m, const Vector<Point>& der1,
                                 const Vector<Point>& der2) {
  const Real dx1 = der1[m][0];
  const Real dx2 = der2[m][0];
  const Real dy1 = der1[m][1];
  const Real dy2 = der2[m][1];
  Tensor<Real, 2> Localdderdx1 = calLocaldderdx1<Order>(pts, crv, m);
  Tensor<Real, 2> Localdderdx2 = calLocaldderdx2<Order>(pts, crv, m);
  Tensor<Real, 2> res(Vec<int, 2>{2, Order + 2});
  for (int j = m - Order / 2; j <= m + Order / 2; j++) {
    res(0, j - m + Order / 2) = -Localdderdx2(0, j - m + Order / 2) * dy1 -
                                dx2 * Localdderdx1(2, j - m + Order / 2) +
                                Localdderdx2(2, j - m + Order / 2) * dx1 +
                                dy2 * Localdderdx1(0, j - m + Order / 2);
    res(1, j - m + Order / 2) = -Localdderdx2(1, j - m + Order / 2) * dy1 -
                                dx2 * Localdderdx1(3, j - m + Order / 2) +
                                Localdderdx2(3, j - m + Order / 2) * dx1 +
                                dy2 * Localdderdx1(1, j - m + Order / 2);
  }
  res(0, Order + 1) =
      -Localdderdx2(0, Order + 1) * dy1 + Localdderdx2(2, Order + 1) * dx1;
  res(1, Order + 1) =
      -Localdderdx2(1, Order + 1) * dy1 + Localdderdx2(3, Order + 1) * dx1;
  return res;
}

// template <int Order>
// Tensor<Real,2> calLocaldfdx(const Vector<Point>& pts, const Crv& crv,
//                                 const int i, const Vector<Point>&
//                             der1, const Vector<Point>& der2, const
//                             Vector<Real>& d2kappa)
// {

// }

template Vector<Real> calArcLength<2>(const Crv& crv);
template Vector<Real> calArcLength<3>(const Crv& crv);
template Vector<Real> calArcLength<4>(const Crv& crv);

// template Vector<Real> calLocalFD2coes<2>(const Vector<Real>&
// arclength, const int i);
// template Vector<Real> calLocalFD2coes<4>(const Vector<Real>&
// arclength, const int i);
// template Vector<Real> calLocalFD2coes<6>(const Vector<Real>&
// arclength, const int i);

template Vector<Real> calLocalFDcoes<2>(const Vector<Real>& arclength,
                                        const int i, const int M);
template Vector<Real> calLocalFDcoes<4>(const Vector<Real>& arclength,
                                        const int i, const int M);
template Vector<Real> calLocalFDcoes<6>(const Vector<Real>& arclength,
                                        const int i, const int M);

template Tensor<Real, 2> calFDcoes<2>(const Vector<Real>& arclength,
                                      const int k);
template Tensor<Real, 2> calFDcoes<4>(const Vector<Real>& arclength,
                                      const int k);
template Tensor<Real, 2> calFDcoes<6>(const Vector<Real>& arclength,
                                      const int k);

// template Tensor<Real,2> calFD1coes<2>(const Vector<Real>& arclength);
// template Tensor<Real,2> calFD1coes<4>(const Vector<Real>& arclength);
// template Tensor<Real,2> calFD1coes<6>(const Vector<Real>& arclength);

// template Tensor<Real,2> calFD2coes<2>(const Vector<Real>& arclength);
// template Tensor<Real,2> calFD2coes<4>(const Vector<Real>& arclength);
// template Tensor<Real,2> calFD2coes<6>(const Vector<Real>& arclength);

template Vector<Point> calDer<2>(const Vector<Point>& datas, const Crv& crv,
                                 const int k);
template Vector<Point> calDer<4>(const Vector<Point>& datas, const Crv& crv,
                                 const int k);
template Vector<Point> calDer<6>(const Vector<Point>& datas, const Crv& crv,
                                 const int k);

template Vector<Real> calDer<2>(const Vector<Real>& datas, const Crv& crv,
                                const int k);
template Vector<Real> calDer<4>(const Vector<Real>& datas, const Crv& crv,
                                const int k);
template Vector<Real> calDer<6>(const Vector<Real>& datas, const Crv& crv,
                                const int k);

template Tensor<Real, 2> calLocaldads1<2>(const Vector<Point>& pts,
                                          const Crv& crv, const int i,
                                          const double dx);
template Tensor<Real, 2> calLocaldads1<4>(const Vector<Point>& pts,
                                          const Crv& crv, const int i,
                                          const double dx);
template Tensor<Real, 2> calLocaldads1<6>(const Vector<Point>& pts,
                                          const Crv& crv, const int i,
                                          const double dx);

template Tensor<Real, 2> calLocaldads2<2>(const Vector<Point>& pts,
                                          const Crv& crv, const int i,
                                          const double dx);
template Tensor<Real, 2> calLocaldads2<4>(const Vector<Point>& pts,
                                          const Crv& crv, const int i,
                                          const double dx);
template Tensor<Real, 2> calLocaldads2<6>(const Vector<Point>& pts,
                                          const Crv& crv, const int i,
                                          const double dx);

template Tensor<Real, 3> caldads1<2>(const Vector<Point>& pts, const Crv& crv,
                                     const double dx);
template Tensor<Real, 3> caldads1<4>(const Vector<Point>& pts, const Crv& crv,
                                     const double dx);
template Tensor<Real, 3> caldads1<6>(const Vector<Point>& pts, const Crv& crv,
                                     const double dx);

template Tensor<Real, 3> caldads2<2>(const Vector<Point>& pts, const Crv& crv,
                                     const double dx);
template Tensor<Real, 3> caldads2<4>(const Vector<Point>& pts, const Crv& crv,
                                     const double dx);
template Tensor<Real, 3> caldads2<6>(const Vector<Point>& pts, const Crv& crv,
                                     const double dx);

// template Vector<Real> calLocaldderdx2<2>(const Vector<Point>& pts, const Crv&
// crv,
//                                          const int i, const int di);
// template Vector<Real> calLocaldderdx2<4>(const Vector<Point>& pts, const Crv&
// crv,
//                                          const int i, const int di);
// template Vector<Real> calLocaldderdx2<6>(const Vector<Point>& pts, const Crv&
// crv,
//                                          const int i, const int di);

template Tensor<Real, 2> calLocaldderdx2<2>(const Vector<Point>& pts,
                                            const Crv& crv, const int i);
template Tensor<Real, 2> calLocaldderdx2<4>(const Vector<Point>& pts,
                                            const Crv& crv, const int i);
template Tensor<Real, 2> calLocaldderdx2<6>(const Vector<Point>& pts,
                                            const Crv& crv, const int i);

template Tensor<Real, 2> calLocaldderdx1<2>(const Vector<Point>& pts,
                                            const Crv& crv, const int i);
template Tensor<Real, 2> calLocaldderdx1<4>(const Vector<Point>& pts,
                                            const Crv& crv, const int i);
template Tensor<Real, 2> calLocaldderdx1<6>(const Vector<Point>& pts,
                                            const Crv& crv, const int i);

template Tensor<Real, 2> calLocaldkappadx<2>(const Vector<Point>& pts,
                                             const Crv& crv, const int m,
                                             const Vector<Point>& der1,
                                             const Vector<Point>& der2);
template Tensor<Real, 2> calLocaldkappadx<4>(const Vector<Point>& pts,
                                             const Crv& crv, const int m,
                                             const Vector<Point>& der1,
                                             const Vector<Point>& der2);
template Tensor<Real, 2> calLocaldkappadx<6>(const Vector<Point>& pts,
                                             const Crv& crv, const int m,
                                             const Vector<Point>& der1,
                                             const Vector<Point>& der2);
