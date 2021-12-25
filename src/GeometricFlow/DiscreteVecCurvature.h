#ifndef _DISCRETEVECCURVATURE_H_
#define _DISCRETEVECCURVATURE_H_

#include "DiscreteVecCoDimOne.h"
#include "FlowHelper.h"

template<int Dim, int Order>
class DiscreteVecCurvature: public DiscreteVecCoDimOne<Dim>
{
private:
  template<class T>
  using Vector = std::vector<T>;

  using Point = Vec<Real, Dim>;
  
public:
  const Vector<Point> operator()(const Vector<Point>& pts, Real t = 0) const;
  
  const Tensor<Real,2> getJacobi(const Vector<Point>& pts, Real t = 0) const;
};

template<int Order>
class DiscreteVecCurvature<2,Order>: public DiscreteVecCoDimOne<2>
{
private:
  template<class T>
  using Vector = std::vector<T>;

  using Point = Vec<Real, 2>;
  
public:
  const Vector<Point> operator()(const Vector<Point>& pts, Real t = 0) const;
  
  const Tensor<Real,2> getJacobi(const Vector<Point>& pts, Real t = 0) const;
};



template<int Order>
const Vector<Point>
DiscreteVecCurvature<2,Order>::operator()(const Vector<Point>& pts, Real t) const
{
  const Curve<2,4> crv = fitCurve<4>(pts,true);
  Vector<Point> res = calDer<Order>(pts,crv,2);
  return res;
}

template<int Order>
const Tensor<Real,2>
DiscreteVecCurvature<2,Order>::getJacobi(const Vector<Point>& pts, Real t) const
{
  const int num = pts.size()-1;
  const Curve<2,4> crv = fitCurve<4>(pts,true);
  Tensor<Real,2> res(Vec<int,2>{num*2,num*2});
  for (int l = 0 ; l < num*2 ; l++)
    for (int m = 0 ; m < num*2 ; m++)
      res(l,m) = 0.0;
  for (int i = 0 ; i < num; i++){
    Tensor<Real,2> Localdderdx2 = calLocaldderdx2<Order>(pts,crv,i);
    for (int j = 0 ; j < Order + 2 ; j++){
      int tmp = i - Order/2 + j;
      if (tmp >= 0 && tmp < num );
      else if (tmp < 0)
        tmp = tmp + num;
      else
        tmp = tmp - num;
      res(i,tmp) = Localdderdx2(0,j);
      res(i,num+tmp) = Localdderdx2(1,j);
      res(num+i,tmp) = Localdderdx2(2,j);
      res(num+i,num+tmp) = Localdderdx2(3,j);
    }
  }
  return res; 
}



#endif

