#ifndef SEGMENTEDREALIZABLESPADJOR_H
#define SEGMENTEDREALIZABLESPADJOR_H

#include <vector>

#include "Core/Curve.h"
#include "Segment.h"
#include "YinSet/OrientedJordanCurve.h"

template <int Order>
class SegmentedRealizableSpadjor;

template <int Dim, int Order>
class YinSet;
template <int Order>
struct CutCellHelper;

template <int Order>
SegmentedRealizableSpadjor<Order> meet(
    const SegmentedRealizableSpadjor<Order>& lhs,
    const SegmentedRealizableSpadjor<Order>& rhs, Real tol);

//============================================================

template <int Order>
class SegmentedRealizableSpadjor {
 public:
  enum { Dim = 2 };
  using rVec = Vec<Real, Dim>;
  using SRS = SegmentedRealizableSpadjor<Order>;
  template <class T>
  using vector = std::vector<T>;

  /// Set up a segmented Spadjor from a Spadjor.
  /// \param aSpadjor
  /// \param tolForSegmentation The tolerance for segmenting the Jordan curves
  /// at the improper intersections. Set to 0.0 if no segmentation is needed.
  SegmentedRealizableSpadjor(
      const vector<OrientedJordanCurve<Dim, Order>>& aSpadjor,
      Real tolForSegmentation = 0.0);
  SegmentedRealizableSpadjor(const vector<Curve<Dim, Order>>& aSpadjor,
                             Real tolForSegmentation = 0.0);
  SegmentedRealizableSpadjor(vector<OrientedJordanCurve<Dim, Order>>&& aSpadjor)
      : orientedJordanCurves(std::move(aSpadjor)){};

  /// Load a segmented Spadjor from the stream. See YinSet::dump().
  /// \param is
  /// \param tolForSegmentation The tolerance for segmenting the Jordan curves
  /// at the improper intersections. Set to 0.0 if no segmentation is needed.
  SegmentedRealizableSpadjor(std::istream& is, Real tolForSegmentation = 0.0);

  // Copy constructor/assignment and move constructor/assignment will go by
  // default.

 public:
  /// Return if the underlying Yin set is bounded.
  bool isBounded(Real tol) const;

  /// Translate the Yin set.
  SRS translate(const rVec& delta) const;

  /// Return the complemententation.
  SRS complement(Real tolForSegmentation = 1e-10) const;

  /// Calculate the meet of two segmented Spadjors.
  friend SRS meet<Order>(const SRS& lhs, const SRS& rhs, Real tol);

  // auto& getOrientedJordanCurves() { return orientedJordanCurves; }
  const auto& getOrientedJordanCurvesRef() const {
    return orientedJordanCurves;
  }

  Real area() const {
    Real volume = 0;
    for (auto& crv : orientedJordanCurves) {
      volume += ::area(crv);
    }
    return volume;
  }

  bool empty() const { return orientedJordanCurves.empty(); }

 protected:
  SegmentedRealizableSpadjor() = default;

  // past Curves and get OrientedJordanCurves
  static void pastSegmentation(
      vector<OrientedJordanCurve<Dim, Order>>& aSpadjor,
      const vector<Curve<Dim, Order>>& crv, Real tolForSegmentation);

  // Cut OrientedJordanCurves and get Curves
  static void applySegmentation(
      const vector<OrientedJordanCurve<Dim, Order>>& aSpadjor,
      vector<Curve<Dim, Order>>& Curves, Real tolForSegmentation);

 protected:
  friend class YinSet<2, Order>;

  vector<OrientedJordanCurve<Dim, Order>> orientedJordanCurves;
};

#endif  // SEGMENTEDREALIZABLESPADJOR_H
