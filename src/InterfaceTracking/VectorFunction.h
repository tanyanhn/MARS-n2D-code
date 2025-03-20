#ifndef _VECTORFUNCTION_H_
#define _VECTORFUNCTION_H_

// #include <Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <vector>

#include "Core/Config.h"
#include "Core/HighPrecisionNumber.h"
#include "Core/Tensor.h"
#include "Core/Vec.h"
#include "Recorder/Timer.h"

template <int Dim>
class VectorFunction {
 private:
  using Point = Vec<Real, Dim>;

 protected:
  using Array = Eigen::Array<HighPrecisionNumber, -1, -1>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  virtual ~VectorFunction(){};

  virtual const Point operator()(const Point& pt, Real t) const = 0;

  virtual const Vector<Real> getJacobi(const Point& pt, Real t) const {
    return Vector<Real>(Dim * Dim);
  }

  virtual const Vector<Point> operator()(const Vector<Point>& pts,
                                         Real t = 0) const {
    Timer timer("VectorFunction::operator()");
    int num = pts.size();
    Vector<Point> vel(num);
    for (int i = 0; i < num; i++) {
      vel[i] = this->operator()(pts[i], t);
    }
    return vel;
  }

  virtual void operator()(const Array& pts, Real t, Array& output) const {}

  static std::string staticClassName() { return "VectorFunction"; }

  /*
  virtual const Tensor<Real, 2> getJacobi(const Vector<Point> &pts, Real t = 0)
  const
  {
      int num = pts.size() - 1;
      Tensor<Real, 2> jacobi(Dim * num);
      jacobi = 0.0;
      Vector<Real> ptjac;

      for (int i = 0; i < num; i++)
      {
          ptjac = getJacobi(pts[i], t);
          for (int j = 0; j < Dim; j++)
          {
              for (int k = 0; k < Dim; k++)
              {
                  jacobi(i + j * num, i + k * num) = ptjac[j * Dim + k];
              }
          }
      }
      return jacobi;
  }
  */

  /*
   virtual const SpMat getJacobi(const Vector<Point> &pts, Real t = 0) const
   {
       int num = pts.size() - 1;
       Vector<Tri> coef;
       coef.reserve(Dim * Dim * num);
       Vector<Real> ptjac;

       for (int i = 0; i < num; i++)
       {
           ptjac = getJacobi(pts[i], t);
           for (int j = 0; j < Dim; j++)
           {
               for (int k = 0; k < Dim; k++)
               {
                   coef.push_back(Tri(i + j * num, i + k * num, ptjac[j * Dim +
   k]));
               }
           }
       }

       SpMat jacobi(num * Dim, num * Dim);
       jacobi.setFromTriplets(coef.begin(), coef.end());
       return jacobi;
   }
   */
};

#endif