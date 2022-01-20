#ifndef YINSET_H
#define YINSET_H

#include "SegmentedRealizableSpadjor.h"

template <int Order>
YinSet<2, Order> intersect(const YinSet<2, Order>& lhs,
                           const YinSet<2, Order>& rhs,
                           Real tol);

template <int Dim, int Order>
struct CurveFactory;
template <int Order>
class YinSet<2, Order> : public SegmentedRealizableSpadjor<Order> {
 public:
  enum { Dim = 2 };
  using SRS = SegmentedRealizableSpadjor<Order>;

  /// Initialize a Yin set from the stream.
  ///
  /// The Hasse diagram is caclulated, but no re-segmentation is applied.
  /// \param tol The tolerance for building the Hasse diagram.
  YinSet(std::istream& is, Real tol);

  /// Initialize a Yin set from a segmented realizable spadjor.
  ///
  /// The pasting map is applied on the collection of segmented curves.
  YinSet(const SRS& segmentedSpadjor, Real tol);

  /// Check if two Yin sets are equal.
  bool equal(const YinSet<Dim, Order>& rhs, Real tol) const;

  /// A node in the Hasse diagram.
  struct Node {
    int depth;  // even number for positive orientation, odd number of
                // negative orientation
    int parent;
    std::vector<int> children;
  };

  /// The Hasse diagram. The last node is the root of the forest.
  std::vector<Node> diagram;

 public:
  /// Get the boundary Jordan curves.
  const std::vector<Curve<Dim, Order>> getBoundaryCycles() const {
    return segmentedCurves;
  }

  /// Return if the Yin set is bounded, based on the Hasse diagram.
  bool isBounded() const { return diagram.back().depth % 2 == 1; }

  /// Get the orientation of the k-th Jordan curve, based on the Hasse
  /// diagram.
  int getOrientation(int k) const {
    return (diagram[k].depth % 2 == 0) ? (1) : (-1);
  }

  /// Return the Betti numbers.
  int getBettiNumber(int rank) const { return bettiNumbers[rank]; }

  /// Return a tabular version of the Hasse diagram.
  std::string getHasseString() const;

  /// Save the Yin set.
  void dump(std::ostream& os) const;

 public:
  /// Calculate the complementation.
  YinSet<2, Order> complement(Real tol) const;

  /// Calculate the intersection of two Yin sets.
  friend YinSet<2, Order> intersect<Order>(const YinSet<2, Order>& lhs,
                                           const YinSet<2, Order>& rhs,
                                           Real tol);

  // YinSet Factory
  friend struct CurveFactory<2, Order>;

 protected:
  ///
  void buildHasse(Real tol);

  int bettiNumbers[2];
  using SRS::segmentedCurves;
};

#endif  // YINSET_H
