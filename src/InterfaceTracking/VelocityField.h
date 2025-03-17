#ifndef _VELOCITYFIELD_H_
#define _VELOCITYFIELD_H_

#include <cmath>

#include "InterfaceTracking/TimeIntegrator.h"
#include "InterfaceTracking/VectorFactory.h"

class CircleShrink : public VectorFunction<2> {
 private:
  using Point = Vec<Real, 2>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  const Point operator()(const Point &pt, Real t) const {
    Vec<Real, 2> y;
    Real l = norm(pt, 2);
    y[0] = -pt[0] / l * v;
    y[1] = -pt[1] / l * v;
    return y;
  }

  const Vector<Real> getJacobi(const Point &pt, Real t) const {
    Vector<Real> ptjac(4);
    Real l = norm(pt, 2);
    ptjac[0] = -v / l;
    ptjac[1] = 0;
    ptjac[2] = 0;
    ptjac[3] = -v / l;
    return ptjac;
  }

 public:
  CircleShrink(Real _v) : v(_v) {}

 private:
  Real v;
};

class SquareShrink : public VectorFunction<2> {
 private:
  using Point = Vec<Real, 2>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  const Point operator()(const Point &pt, Real t) const {
    Real l = std::max(std::abs(pt[0]), std::abs(pt[1]));
    Point y;
    y[0] = pt[0] / l;
    y[1] = pt[1] / l;
    return y;
  }

  const Vector<Real> getJacobi(const Point &pt, Real t) const {
    Vector<Real> ptjac(4);
    return ptjac;
  }
};

class Translation : public VectorFunction<2> {
 private:
  using Point = Vec<Real, 2>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  const Point operator()(const Point &pt, Real t) const {
    Vec<Real, 2> y;
    y[0] = v1;
    y[1] = v2;
    return y;
  }

  const Vector<Real> getJacobi(const Point &pt, Real t) const {
    Vector<Real> ptjac(4);
    ptjac[0] = 0;
    ptjac[1] = 0;
    ptjac[2] = 0;
    ptjac[3] = 0;
    return ptjac;
  }

 public:
  Translation(Real u1, Real u2) : v1(u1), v2(u2){};

 private:
  Real v1, v2;
};

class Rotation : public VectorFunction<2> {
 private:
  using Point = Vec<Real, 2>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  const Point operator()(const Point &pt, Real t) const {
    Vec<Real, 2> y;
    y[0] = om * (rc1 - pt[1]);
    y[1] = -om * (rc2 - pt[0]);
    return y;
  }

  const Vector<Real> getJacobi(const Point &pt, Real t) const {
    Vector<Real> ptjac(4);
    ptjac[0] = 0;
    ptjac[1] = -om;
    ptjac[2] = om;
    ptjac[3] = 0;
    return ptjac;
  }

 public:
  Rotation(Real c1, Real c2, Real _om) : rc1(c1), rc2(c2), om(_om){};

 private:
  Real rc1, rc2, om;
};

class RevRotation : public VectorFunction<2> {
 private:
  using Point = Vec<Real, 2>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  const Point operator()(const Point &pt, Real t) const {
    Point y;
    y[0] = om * (rc2 - pt[1]) * cos(M_PI * t / T);
    y[1] = -om * (rc1 - pt[0]) * cos(M_PI * t / T);
    return y;
  }

  const Vector<Real> getJacobi(const Point &pt, Real t) const {
    Vector<Real> ptjac(4);
    ptjac[0] = 0;
    ptjac[1] = -om * cos(M_PI * t / T);
    ptjac[2] = om * cos(M_PI * t / T);
    ptjac[3] = 0;
    return ptjac;
  }

 public:
  RevRotation(Real c1, Real c2, Real _om, Real _T)
      : rc1(c1), rc2(c2), om(_om), T(_T){};

 private:
  Real rc1, rc2, om, T;
};

class Vortex : public VectorFunction<2> {
 private:
  using Point = Vec<Real, 2>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  const Point operator()(const Point &pt, Real t) const override {
    Point y;
    y[0] =
        cos(M_PI * t / T) * pow(sin(M_PI * pt[0]), 2) * sin(2 * M_PI * pt[1]);
    y[1] =
        -cos(M_PI * t / T) * pow(sin(M_PI * pt[1]), 2) * sin(2 * M_PI * pt[0]);
    return y;
  }

  void operator()(const Eigen::ArrayXXd &pts, Real t,
                  Eigen::ArrayXXd &output) const override {
    Timer timer("vectorFunction::operator()");
    using namespace Eigen;
    // const long N = pts.cols();

    const double scale = cos(M_PI * t / T);
    // 预先计算公共因子x = M_PI * pts.array()
    const ArrayXXd x = M_PI * pts.array();
    const ArrayXXd cos_term = x.cos();
    const ArrayXXd sin_term = x.sin();
    const ArrayXXd squared_sin = sin_term.square();
    const ArrayXXd double_sin = 2 * sin_term * cos_term;
    // const ArrayXXd double_sin = (2 * x).sin();

    // ArrayXXd result(2, N);
    output.row(0) = scale * squared_sin.row(0) * double_sin.row(1);
    output.row(1) = -scale * squared_sin.row(1) * double_sin.row(0);
  }

  const Vector<Real> getJacobi(const Point &pt, Real t) const override {
    Vector<Real> ptjac(4);
    ptjac[0] = cos(M_PI * t / T) * 2 * M_PI * sin(M_PI * pt[0]) *
               cos(M_PI * pt[0]) * sin(2 * M_PI * pt[1]);
    ptjac[1] = cos(M_PI * t / T) * 2 * M_PI * pow(sin(M_PI * pt[0]), 2) *
               cos(2 * M_PI * pt[1]);
    ptjac[2] = -cos(M_PI * t / T) * 2 * M_PI * sin(M_PI * pt[1]) *
               cos(M_PI * pt[1]) * sin(2 * M_PI * pt[0]);
    ptjac[3] = -cos(M_PI * t / T) * 2 * M_PI * pow(sin(M_PI * pt[1]), 2) *
               cos(2 * M_PI * pt[0]);
    return ptjac;
  }

 public:
  Vortex(Real _T) : T(_T){};
  static std::string staticClassName() { return "Vortex"; }

 private:
  Real T;
};
REGISTER_CLASS(Vortex, Real) // 注册参数类型

class Deformation : public VectorFunction<2> {
 private:
  using Point = Vec<Real, 2>;

  template <class T>
  using Vector = std::vector<T>;

 public:
  const Point operator()(const Point &pt, Real t) const override{
    Point y;
    y[0] = cos(M_PI * t / T) * sin(n * M_PI * (pt[0] + 0.5)) *
           sin(n * M_PI * (pt[1] + 0.5));
    y[1] = cos(M_PI * t / T) * cos(n * M_PI * (pt[0] + 0.5)) *
           cos(n * M_PI * (pt[1] + 0.5));
    return y;
  }

  void operator()(const Eigen::ArrayXXd &pts, Real t,
                  Eigen::ArrayXXd &output) const override {
    Timer timer("vectorFunction::operator()");
    using namespace Eigen;

    const double scale = cos(M_PI * t / T);
    const double n_pi = n * M_PI;

    // 计算变换后的坐标
    ArrayXd x = n_pi * (pts.row(0).array() + 0.5);
    ArrayXd y = n_pi * (pts.row(1).array() + 0.5);

    // 计算sin和cos值
    ArrayXd sin_x = x.sin();
    ArrayXd sin_y = y.sin();
    ArrayXd cos_x = x.cos();
    ArrayXd cos_y = y.cos();

    // 组合结果
    output.row(0) = scale * (sin_x * sin_y);
    output.row(1) = scale * (cos_x * cos_y);
  }

  const Vector<Real> getJacobi(const Point &pt, Real t) const override{
    Vector<Real> ptjac(4);
    ptjac[0] = cos(M_PI * t / T) * n * M_PI * cos(n * M_PI * (pt[0] + 0.5)) *
               sin(n * M_PI * (pt[1] + 0.5));
    ptjac[1] = cos(M_PI * t / T) * n * M_PI * sin(n * M_PI * (pt[0] + 0.5)) *
               cos(n * M_PI * (pt[1] + 0.5));
    ptjac[2] = -cos(M_PI * t / T) * n * M_PI * cos(n * M_PI * (pt[0] + 0.5)) *
               sin(n * M_PI * (pt[1] + 0.5));
    ptjac[3] = -cos(M_PI * t / T) * n * M_PI * sin(n * M_PI * (pt[0] + 0.5)) *
               cos(n * M_PI * (pt[1] + 0.5));
    return ptjac;
  }

 public:
  Deformation(Real _T, int _n = 4) : T(_T), n(_n){};
  static std::string staticClassName() { return "Deformation"; }

 private:
  Real T;
  int n;
};
REGISTER_CLASS(Deformation, Real, int) // 注册参数类型

#endif