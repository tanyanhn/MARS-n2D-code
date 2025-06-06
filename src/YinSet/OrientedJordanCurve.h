#ifndef ORIENTEDJORDANCURVE_H
#define ORIENTEDJORDANCURVE_H

#include <memory>
#include <string>
#include <vector>

#include "Core/Curve.h"
#include "Core/VecCompare.h"
#include "SimplicialComplex.h"

template <int Dim, int Order>
class YinSet;

template <int Dim, int Order>
class OrientedJordanCurve : public Curve<Dim, Order> {
 public:
  // (i, j), i-th Curve, j-th knot point.
  using PointIndex = std::pair<size_t, size_t>;
  using Vertex = Simplex<PointIndex>::Vertex;
  // constructors
  OrientedJordanCurve() = default;

  OrientedJordanCurve(const Curve<Dim, Order>& curve)
      : Curve<Dim, Order>(curve) {
    // makeSelfMonotonic(distTol());
    VecCompare<Real, Dim> vCmp(distTol());
    assert(curve.empty() || vCmp(curve.startpoint(), curve.endpoint()) == 0 &&
                                "Initial OrientedJordanCurve must maintain "
                                "startpoint == endpoint.");
  }
  virtual ~OrientedJordanCurve() = default;

  // virtual constructor
  virtual void define(const std::string& parameters);

  virtual void define(const std::string& parameters,
                      SimplicialComplex<Vertex>& kinks);

  virtual void define(std::istream& is, SimplicialComplex<Vertex>& kinks);

  // initial function
  void define(const std::vector<Vec<Real, Dim>>& points,
              const SimplicialComplex<Vertex>& kinks);
  void define(const std::vector<Vec<Real, Dim>>& points,
              const vector<Vertex>& indexes);

  void makeSelfMonotonic(Real tol) { *this = this->makeMonotonic(tol).first; }

  // locate the point in the OrientedJordanCurve.
  enum class LocateResult {
    // point is inside the curve.
    INSIDE,
    // point is outside the curve.
    OUTSIDE,
    // point is on the curve.
    ON_CURVE
  };
  auto locate(const Vec<Real, Dim>& point, Real tol = distTol()) const
      -> LocateResult;

  template <int, int>
  friend class YinSet;
};

//==========================================================

template <int Order>
struct Circle : public OrientedJordanCurve<2, Order> {
  using Vertex = typename OrientedJordanCurve<2, Order>::Vertex;

  // constructors
  Circle() = default;

  // virtual constructor
  void define(const std::string& parameters);

  void define(const std::string& parameters, SimplicialComplex<Vertex>& kinks);

  void define(std::istream& is, SimplicialComplex<Vertex>& kinks);
};

template <int Order>
struct Rectangle : public OrientedJordanCurve<2, Order> {
  using Vertex = typename OrientedJordanCurve<2, Order>::Vertex;

  // constructors
  Rectangle() = default;

  // virtual constructor
  void define(const std::string& parameters);

  void define(const std::string& parameters, SimplicialComplex<Vertex>& kinks);

  void define(std::istream& is, SimplicialComplex<Vertex>& kinks);
};

//=========================================================================================

template <int Dim, int Order>
struct CurveFactory {
  CurveFactory() = default;
  std::unique_ptr<OrientedJordanCurve<Dim, Order>> createCurve(
      const std::string& parameters);
};

template <int Order>
struct CurveFactory<2, Order> {
  using PointIndex = typename YinSet<2, Order>::PointIndex;
  using Vertex = typename YinSet<2, Order>::Vertex;

  CurveFactory() = default;

  // factory function
  std::unique_ptr<OrientedJordanCurve<2, Order>> createCurve(
      const std::string& parameters);

  std::unique_ptr<OrientedJordanCurve<2, Order>> createCurve(
      const std::string& parameters, SimplicialComplex<Vertex>& kinks);

  std::unique_ptr<OrientedJordanCurve<2, Order>> createCurve(
      std::istream& is, SimplicialComplex<Vertex>& kinks);

  YinSet<2, Order> createYinSet(const std::vector<std::string>& parameters);
};

#endif  // !ORIENTEDJORDANCURVE_H
