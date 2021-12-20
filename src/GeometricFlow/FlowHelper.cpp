#include "FlowHelper.h"
#include "Core/numlib.h"
#include <cmath>

template <int Order>
Real calArcLength(const Polynomial<4, Point> &poly, Real start, Real end)
{
    auto dxdt = getComp(poly, 0).der();
    auto dydt = getComp(poly, 1).der();
    auto ndsdt = [&](Real t)
    {
        Vec<Real, 2> ds{dxdt(t), dydt(t)};
        return norm(ds, 2);
    };
    return quad<Order>(ndsdt, start, end);
}


template <int Order>
Vector<Real> calArcLength(const Crv& crv)
{
  const Vector<Polynomial<4, Point> >& polys = crv.getPolys();
  const Vector<Real>& knots = crv.getKnots();
  const int num = knots.size();
  Vector<Real> res(num);
  res[0] = 0.0;
  for (int i = 1 ; i < num ; i++)
    res[i] = res[i-1] + calArcLength<Order>(polys[i-1],0,knots[i]-knots[i-1]);
  return res;
}


template <int Order>
Tensor<Real,2> calFD1coes(const Vector<Real>& arclength)
{
  assert(Order%2 == 0);
  const int num = arclength.size();
  Tensor<Real,2> res(Vec<int,2> {Order+1,num-1});
  for (int i = 0 ; i < num - 1 ; i++){
    for (int j = i - Order/2 ; j <= i + Order/2 ; j++){
      Real tmpj = 0;
      if (j >= 0 && j < num)
        tmpj = arclength[j];
      else if (j < 0)
        tmpj = arclength[j+num-1] - arclength[num-1];
      else
        tmpj = arclength[num-1] + arclength[j-num+1];
      Polynomial<Order+1,Real> poly;
      poly[0] = 1.0;
      for (int k = i - Order/2 ; k <= i + Order/2 ; k++){
        if (k != j){
          Real tmpk = 0;
          if (k >= 0 && k < num)
            tmpk = arclength[k];
          else if (k < 0)
            tmpk = arclength[k+num-1] - arclength[num-1];
          else
            tmpk = arclength[num-1] + arclength[k-num+1];
          Polynomial<2,Real>
            tmppoly{-(tmpk-arclength[i])/(tmpj-tmpk),1/(tmpj-tmpk)};
          auto polyl = poly*tmppoly;
          for (int m = 0 ; m < std::min(k - i + Order/2 + 2,Order+1); m++)
            poly[m] = polyl[m];
        }
      }
      auto dpoly = poly.der();
      res(j-i+Order/2,i) = dpoly(0.0);
    }
  }
  return res;
}


template <int Order>
Tensor<Real,2> calFD2coes(const Vector<Real>& arclength)
{
  assert(Order%2 == 0);
  const int num = arclength.size();
  Tensor<Real,2> res(Vec<int,2> {Order+2,num-1});
  for (int i = 0 ; i < num - 1 ; i++){
    for (int j = i - Order/2 ; j <= i + Order/2 + 1 ; j++){
      Real tmpj = 0;
      if (j >= 0 && j < num)
        tmpj = arclength[j];
      else if (j < 0)
        tmpj = arclength[j+num-1] - arclength[num-1];
      else
        tmpj = arclength[num-1] + arclength[j-num+1];
      Polynomial<Order+2,Real> poly;
      poly[0] = 1.0;
      for (int k = i - Order/2 ; k <= i + Order/2 + 1 ; k++){
        if (k != j){
          Real tmpk = 0;
          if (k >= 0 && k < num)
            tmpk = arclength[k];
          else if (k < 0)
            tmpk = arclength[k+num-1] - arclength[num-1];
          else
            tmpk = arclength[num-1] + arclength[k-num+1];
          Polynomial<2,Real>
            tmppoly{-(tmpk-arclength[i])/(tmpj-tmpk),1/(tmpj-tmpk)};
          auto polyl = poly*tmppoly;
          for (int m = 0 ; m < std::min(k - i + Order/2 + 2,Order+2); m++)
            poly[m] = polyl[m];
        }
      }
      auto ddpoly = poly.der().der();
      res(j-i+Order/2,i) = ddpoly(0.0);
    }
  }
  return res;
}

template <int Order, class T>
Vector<T> calDer(const Vector<T>& datas, const Crv& crv, const int k)
{
  assert(Order%2 == 0);
  Vector<Real> arclength = calArcLength<Order/2+1>(crv);
  Tensor<Real,2> FDcoes;
  switch (k){
  case 1:
    FDcoes = calFD1coes<Order>(arclength);
    break;
  case 2:
    FDcoes = calFD2coes<Order>(arclength);
    break;
  default:
    throw Vector<T>();
  }
  const int num = datas.size();
  Vector<T> res(num);
  const int numterms = (FDcoes.size())[0];
  for (int i = 0 ; i < num - 1 ; i++){
    for (int j = 0 ; j < numterms ; j++){
      int tmp = 0;
      const int index = i - Order/2 + j;
      if (index >= 0 && index < num)
        tmp = index;
      else if (index < 0)
        tmp = index + num - 1;
      else
        tmp = index - num + 1;
      res[i] = res[i] + datas[tmp]*FDcoes(j,i);
    }
  }
  res[num-1] = res[0];
  return res;
}




template Vector<Real> calArcLength<2>(const Crv& crv);
template Vector<Real> calArcLength<3>(const Crv& crv);
template Vector<Real> calArcLength<4>(const Crv& crv);

template Tensor<Real,2> calFD1coes<2>(const Vector<Real>& arclength);
template Tensor<Real,2> calFD1coes<4>(const Vector<Real>& arclength);
template Tensor<Real,2> calFD1coes<6>(const Vector<Real>& arclength);

template Tensor<Real,2> calFD2coes<2>(const Vector<Real>& arclength);
template Tensor<Real,2> calFD2coes<4>(const Vector<Real>& arclength);
template Tensor<Real,2> calFD2coes<6>(const Vector<Real>& arclength);

template Vector<Point> calDer<2>(const Vector<Point>& datas, const 
Crv& crv, const int k);
template Vector<Point> calDer<4>(const Vector<Point>& datas, const 
Crv& crv, const int k);
template Vector<Point> calDer<6>(const Vector<Point>& datas, const 
Crv& crv, const int k);

template Vector<Real> calDer<2>(const Vector<Real>& datas, const 
Crv& crv, const int k);
template Vector<Real> calDer<4>(const Vector<Real>& datas, const 
Crv& crv, const int k);
template Vector<Real> calDer<6>(const Vector<Real>& datas, const 
Crv& crv, const int k);
