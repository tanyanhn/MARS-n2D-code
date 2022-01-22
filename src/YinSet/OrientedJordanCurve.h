#ifndef ORIENTEDJORDANCURVE_H
#define ORIENTEDJORDANCURVE_H

#include <memory>
#include <sstream>
#include <vector>
#include "Core/Curve.h"
#include "SimplicialComplex.h"
#include "YinSet/YinSet.h"

template <int Dim, int Order>
class YinSet;

template <int Dim, int Order>
class OrientedJordanCurve : public Curve<Dim, Order> {
 public:
  // constructors
  OrientedJordanCurve() = default;

  virtual void define(const std::string& parameters);

 protected:
  void define(const std::vector<Vec<Real, Dim>>& points,
              const SimplicialComplex& kinks);
};

template <int Order>
struct Circle : public OrientedJordanCurve<2, Order> {
  // constructors
  Circle() = default;

  void define(const std::string& parameters);
};

template <int Order>
struct Rectangle : public OrientedJordanCurve<2, Order> {
  // constructors
  Rectangle() = default;

  void define(const std::string& parameters);
};

template <int Dim, int Order>
struct CurveFactory {
  CurveFactory() = default;
  std::unique_ptr<OrientedJordanCurve<Dim, Order>> createCurve(
      const std::string& parameters);
  // YinSet<Dim, Order> createYinSet(const std::vector<std::string>&
  // parameters);
};

template <int Order>
struct CurveFactory<2, Order> {
  CurveFactory() = default;
  std::unique_ptr<OrientedJordanCurve<2, Order>> createCurve(
      const std::string& parameters);
  YinSet<2, Order> createYinSet(const std::vector<std::string>& parameters);
};

#endif  // !ORIENTEDJORDANCURVE_H
