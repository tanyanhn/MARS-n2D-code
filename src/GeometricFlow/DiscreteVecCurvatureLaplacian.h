#ifndef _DISCRETEVECCURVATURELAPLACIAN_H_
#define _DISCRETEVECCURVATURELAPLACIAN_H_

#include "DiscreteVecCoDimOne.h"
#include "FlowHelper.h"

template<int Dim, int Order>
class DiscreteVecCurvatureLaplacian: public DiscreteVecCoDimOne<Dim>
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
class DiscreteVecCurvatureLaplacian<2,Order>
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
DiscreteVecCurvatureLaplacian<2,Order>::operator()(const
                                                   Vector<Point>& pts, Real t) const{
  const int num = pts.size();
  const Curve<2,4> crv = fitCurve<4>(pts,true);
  Vector<Point> der1 = calDer<Order+2>(pts,crv,1);
  Vector<Point> der2 = calDer<Order+2>(pts,crv,2);
  Vector<Real> kappa = Vector<Real>(num);
  for (int i = 0 ; i < num ; i++){
    kappa[i] = -der2[i][0]*der1[i][1]+der2[i][1]*der1[i][0];
    //std::cout << kappa[i]-1/0.15 << std::endl;
  }
  Vector<Real> d2kappa = calDer<Order+2>(kappa,crv,2);
  Vector<Point> res = Vector<Point>(num);
  for (int i = 0 ; i < num ; i++){
    res[i][0] = d2kappa[i]*der1[i][1];
    res[i][1] = -d2kappa[i]*der1[i][0];
  }
  return res;
}

template<int Order>
const Tensor<Real,2>
DiscreteVecCurvatureLaplacian<2,Order>::getJacobi(const
                                                  Vector<Point>& pts, Real t) const
{  
  return Tensor<Real,2>();
}



#endif
