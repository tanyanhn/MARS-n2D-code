#ifndef _VECTORONHYPERSURFACE_H_
#define _VECTORONHYPERSURFACE_H_

#include "Core/Config.h"
#include "Core/Vec.h"
#include "Core/Tensor.h"
#include <vector>

template <int Dim>
class VectorOnHypersurface
{
private:
  template<class T>
  using Vector = std::vector<T>;
  
  using Point = Vec<Real, Dim>;
 
public:
  virtual const Vector<Point> operator()(const Vector<Point>& pts,
                                         Real t = 0) const = 0;
  
  virtual const Tensor<Real,2> getJacobi(const Vector<Point>& pts,
                                         Real t = 0) const{return Tensor<Real,2>();}
};


#endif
