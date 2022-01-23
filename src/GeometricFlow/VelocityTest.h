#ifndef _VELOCITYTEST_H_
#define _VELOCITYTEST_H_

#include "FlowHelper.h"
#include "InterfaceTracking/TimeIntegrator.h"
#include <algorithm>
#include <numeric>
#include <cmath>

template<class T>
using Vector = std::vector<T>;

using Point = Vec<Real, 2>;

using Crv = Curve<2, 4>;

// class CircleInnerNormal : public VectorFunction<2>
// {
// private:
//   template<class T>
//   using Vector = std::vector<T>;

//   using Point = Vec<Real, 2>;
  
//   const Point operator()(Point pt, Real t) const
//   {
//     Point y =  (pt - center)*(-1.0/(radio*radio));
//     return y;
//   }
  
//   const Vector<Real> getJacobi(Point pt, Real t) const
//   {
//     Vector<Real> ptjac(4);
//     const double tmp = pow(pt[0]*pt[0]+pt[1]*pt[1],1.5);
//     ptjac[0] = -1.0/radio * (pt[1]*pt[1]/tmp);
//     ptjac[1] = -1.0/radio * (-pt[0]*pt[1]/tmp);
//     ptjac[2] = ptjac[1];
//     ptjac[3] = -1.0/radio * (pt[0]*pt[0]/tmp);
//     return ptjac;
//   }

// public:
//   CircleInnerNormal(Point Center, Real Radio) : center(Center), radio(Radio){};

// private:
//   Point center;
//   Real radio;
// };


Real OneNormVelocity(const Vector<Point>& v1, const Vector<Point>& v2,
const double h)
{
  const int num1 = v1.size();
  const int num2 = v2.size();
  assert(num1 == num2);
  Vector<Real> diff;
  for (int i = 0 ; i < num1 ; i++)
    diff.push_back(norm(v1[i]-v2[i],2));
  Real sum = std::accumulate(diff.begin(), diff.end(), 0.0);
  return sum*h;
}
  
Real MaxNormVelocity(const Vector<Point>& v1, const Vector<Point>& v2)
{
  const int num1 = v1.size();
  const int num2 = v2.size();
  assert(num1 == num2);
  Vector<Real> diff;
  for (int i = 0 ; i < num1 ; i++)
    diff.push_back(norm(v1[i]-v2[i],2));
  Real max =*std::max_element(diff.begin(),diff.end());
  return max;
}

Vector<Point> Extrapolation(const Vector<Point>& v1, const
Vector<Point>& v2)
{
  const int num1 = v1.size();
  const int num2 = v2.size();
  assert((num2-1)%(num1-1) == 0);
  const int r = (num2-1)/(num1-1);
  Vector<Point> ref;
  for (int i = 0 ; i < num1 - 1 ; i++)
    ref.push_back(v2[i*r]);
  ref.push_back(v2[0]);
  return ref;
}

Vector<Point> CircleCurvature(const Real radio, const int n)
{
  Vector<Point> pts;
  pts.push_back({-1.0/radio, 0.0});
  for (int i = 1; i < n; i++)
    {
      pts.push_back({-1.0/radio * cos(2 * M_PI / n * i), -1.0/radio * sin(2 * M_PI / n * i)});
    }
  pts.push_back({-1.0/radio, 0.0});
  return pts;
}


Vector<Point> CircleCurvatureLaplacian(const Real radio, const int n)
{
  Vector<Point> pts;
  for (int i = 0; i <= n; i++)
    pts.push_back({0.0, 0.0});
  return pts;
}




#endif
