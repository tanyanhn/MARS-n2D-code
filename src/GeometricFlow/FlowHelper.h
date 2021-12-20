#ifndef _FLOWHELPER_H_
#define _FLOWHELPER_H_

#include "Core/Curve.h"
#include "Core/Tensor.h"

template<class T>
using Vector = std::vector<T>;

using Point = Vec<Real, 2>;

using Crv = Curve<2, 4>;

template <int Order>
Vector<Real> calArcLength(const Crv& crv);

template <int Order>
Tensor<Real,2> calFD1coes(const Vector<Real>& arclength);

template <int Order>
Tensor<Real,2> calFD2coes(const Vector<Real>& arclength);

//crv must be the curve fitted by points lying in datas, and k must be
//1 or 2
template <int Order, class T>
Vector<T> calDer(const Vector<T>& datas, const Crv& crv, const int k);
 








#endif
