#ifndef _DISCRETEVECCODIMONE_H_
#define _DISCRETEVECCODIMONE_H_

#include "InterfaceTracking/VectorFunction.h"



template <int Dim>
class DiscreteVecCoDimOne: public VectorFunction<Dim>
{
private:
  template<class T>
  using Vector = std::vector<T>;
  
  using Point = Vec<Real, Dim>;
  
  virtual const Point operator()(Point pt, Real t) const{throw pt;}
  
  virtual const Vector<Real> getJacobi(Point pt, Real t) const {throw pt;}
public:
  virtual const Vector<Point> operator()(const Vector<Point>& pts,
                                         Real t = 0) const = 0;
  
  virtual const Tensor<Real,2> getJacobi(const Vector<Point>& pts,
                                         Real t = 0) const{return
                                         Tensor<Real,2>();}
};


#endif
