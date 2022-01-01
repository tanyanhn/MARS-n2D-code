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
Vector<Real> calLocalFDcoes(const Vector<Real>& arclength, const int
i, const int M);

template <int Order>
Tensor<Real,2> calFDcoes(const Vector<Real>& arclength, const int k);

//crv must be the curve fitted by points lying in datas, and k must be
//1 or 2
template <int Order, class T>
Vector<T> calDer(const Vector<T>& datas, const Crv& crv, const int k);


//The following functions are used to calculate Jacobi matrix.


// Input: $(x_j,y_j)$ and $(x_{j+1},y_{j+1}$
// Output:
//$$(\frac{\partial s_j}{\partial x_j},\frac{\partial s_j}{\partial x_{j+1},
//  \frac{\partial s_j}{\partial y_j},\frac{\partial s_j}{\partial y_{j+1}})$$
Vector<Real> calLocaldsdx(const Point p1, const Point p2);

// For all adjacent points
Tensor<Real,2> caldsdx(const Vector<Point>& pts);


// Input: tracking points, corresponding curve and index i
// Output:
//$$(\frac{\partial \alpha_{i,i-r/2}}{\partial s_{i-r/2}},
//   \frac{\partial \alpha_{i,i-r/2}}{\partial s_{i-r/2+1}},
//   \cdots,
//   \frac{\partial \alpha_{i,i-r/2}}{\partial s_{i+r/2}});
//   \frac{\partial \alpha_{i,i-r/2+1}}{\partial s_{i-r/2}},
//   \cdots,
//   \frac{\partial \alpha_{i,i-r/2+1}}{\partial s_{i+r/2}});
//   \cdots,
//   \frac{\partial \alpha_{i,i+r/2}}{\partial s_{i+r/2}})$$
// Here r = Order, "1" and "2" in function name denote order of derivative.
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

// Input:tracking points, corresponding curve, index i and di.
// Output:
//$$(\frac{\partial f_{s,r}(x_i)}{\partial x_{i+di}},
//   \frac{\partial f_{s,r}(x_i)}{\partial y_{i+di}},
//   \frac{\partial f_{s,r}(y_i)}{\partial x_{i+di}},
//   \frac{\partial f_{s,r}(y_i)}{\partial y_{i+di}})$$
// Here r = Order, "1" and "2" in function name denote order of derivative.
// Notice when di <= -r/2-1 or di >= r/2+2, the output is all zeros.

//template <int Order>
//Vector<Real> calLocaldderdx2(const Vector<Point>& pts, const Crv& crv,
//                            const int i, const int di);

// For all di in [-r/2,r/2+1], output a 4*(r+2) Tensor.
template <int Order>
Tensor<Real,2> calLocaldderdx2(const Vector<Point>& pts, const Crv& crv,
                               const int i);


//The following functions are used to calculate Jacobi matrix of
//surface diffusion

template <int Order>
Tensor<Real,2> calLocaldderdx1(const Vector<Point>& pts, const Crv& crv,
                               const int i);

// Input: tracking points, corresponding curve, index m, derivative
// value at each tracking point.
// Output:
// $$(\frac{\partial \kappa_m}{\partial x_{m-r/2}},
//    \frac{\partial \kappa_m}{\partial x_{m-r/2+1}},
//    \cdots,
//    \frac{\partial \kappa_m}{\partial x_{m+r/2+1}};
//    \frac{\partial \kappa_m}{\partial y_{m-r/2}},
//    \frac{\partial \kappa_m}{\partial y_{m-r/2+1}},
//    \cdots,
//    \frac{\partial \kappa_m}{\partial y_{m+r/2+1}})$$
// Here r = Order, and output is a 2-by-(Order+2) tensor.
template <int Order>
Tensor<Real,2> calLocaldkappadx(const Vector<Point>& pts, const Crv& crv,
                                const int m, const Vector<Point>&
                                der1, const Vector<Point>& der2);


#endif
