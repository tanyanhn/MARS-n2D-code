#ifndef _DISCRETEVECCURVATURELAPLACIANNEW_H_
#define _DISCRETEVECCURVATURELAPLACIANNEW_H_

#include "DiscreteVecCoDimOne.h"
#include "DiscreteVecCurvatureLaplacianOld.h" 
#include "FlowHelper.h"

template<int Dim, int Order>
class DiscreteVecCurvatureLaplacianNew: public DiscreteVecCoDimOne<Dim>
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
class DiscreteVecCurvatureLaplacianNew<2,Order>: public DiscreteVecCoDimOne<2>
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
DiscreteVecCurvatureLaplacianNew<2,Order>::operator()(const
                                                   Vector<Point>& pts, Real t) const{
  const int num = pts.size();
  const Curve<2,4> crv = fitCurve<4>(pts,true);
  Vector<Point> der1 = calDer<Order>(pts,crv,1);
  Vector<Point> der2 = calDer<Order>(pts,crv,2);
  Vector<Point> der3 = calDer<Order>(pts,crv,3);
  Vector<Point> der4 = calDer<Order>(pts,crv,4);
  Vector<Real> tmp(num);
  for (int i = 0 ; i < num ; i++){
    tmp[i] =
    -der4[i][0]*der1[i][1]-2*der3[i][0]*der2[i][1]-der2[i][0]*der3[i][1]
    +
    der4[i][1]*der1[i][0]+2*der3[i][1]*der2[i][0]+der2[i][1]*der3[i][0];
  }
  Vector<Point> res(num);
  for (int i = 0 ; i < num ; i++){
    res[i][0] = tmp[i]*der1[i][1];
    res[i][1] = -tmp[i]*der1[i][0];
  }
  return res;
}

template<int Order>
const Tensor<Real,2>
DiscreteVecCurvatureLaplacianNew<2,Order>::getJacobi(const
                                                  Vector<Point>& pts, Real t) const
{
  DiscreteVecCurvatureLaplacianOld<2,Order> SDF;
  return SDF(pts,t);
}



#endif
