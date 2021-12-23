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

//Input: $(x_j,y_j)$ and $(x_{j+1},y_{j+1}$
//Output:
//$$(\frac{\partial s_j}{\partial x_j},\frac{\partial s_j}{\partial x_{j+1},
//  \frac{\partial s_j}{\partial y_j},\frac{\partial s_j}{\partial y_{j+1}})$$
Vector<Real> calLocaldsdx(const Point p1, const Point p2);

//For all adjacent points
Tensor<Real,2> caldsdx(const Vector<Point>& pts);


//Input: tracking points , corresponding curve and index i
//Output:
//$$(\frac{\partial \alpha_{i,i-r/2}}{\partial s_{i-r/2}},
//   \frac{\partial \alpha_{i,i-r/2}}{\partial s_{i-r/2+1}},
//   \cdots,
//   \frac{\partial \alpha_{i,i-r/2}}{\partial s_{i+r/2}});
//   \frac{\partial \alpha_{i,i-r/2+1}}{\partial s_{i-r/2}},
//   \cdots,
//   \frac{\partial \alpha_{i,i-r/2+1}}{\partial s_{i+r/2}});
//   \cdots,
//   \frac{\partial \alpha_{i,i+r/2}}{\partial s_{i+r/2}})$$
//Here r = Order.
template <int Order>
Tensor<Real,2> calLocaldads1(const Vector<Point>& pts, const Crv& crv,
const int i , const double dx);

template <int Order>
Tensor<Real,2> calLocaldads2(const Vector<Point>& pts, const Crv& crv,
const int i , const double dx);

//For all index
template <int Order>
Tensor<Real,3> caldads1(const Vector<Point>& pts, const Crv& crv,
                        const double dx);

template <int Order>
Tensor<Real,3> caldads2(const Vector<Point>& pts, const Crv& crv,
                        const double dx);

#endif
