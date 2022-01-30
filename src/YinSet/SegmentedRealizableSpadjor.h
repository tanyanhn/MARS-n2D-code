#ifndef SEGMENTEDREALIZABLESPADJOR_H
#define SEGMENTEDREALIZABLESPADJOR_H

#include <string>
#include <vector>
#include "Core/Curve.h"
#include "Segment.h"
#include "YinSet/OrientedJordanCurve.h"

template <int Order>
class SegmentedRealizableSpadjor;

template <int Dim, int Order>
class YinSet;

template <int Order>
SegmentedRealizableSpadjor<Order> meet(
    const SegmentedRealizableSpadjor<Order>& lhs,
    const SegmentedRealizableSpadjor<Order>& rhs,
    Real tol);

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
  SRS complement() const;

  /// Calculate the meet of two segmented Spadjors.
  friend SRS meet<Order>(const SRS& lhs, const SRS& rhs, Real tol);

 protected:
  SegmentedRealizableSpadjor() = default;

  void applySegmentation(
      const vector<OrientedJordanCurve<Dim, Order>>& aSpadjor,
      Real tolForSegmentation = 0.0);

 protected:
  friend class YinSet<2, Order>;

  vector<Curve<Dim, Order>> Curves;
  vector<OrientedJordanCurve<Dim, Order>> segmentedCurves;
};

#endif  // SEGMENTEDREALIZABLESPADJOR_H
