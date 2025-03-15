#ifndef YINSET_H
#define YINSET_H

#include "YinSet/OrientedJordanCurve.h"
#include "YinSet/SegmentedRealizableSpadjor.h"
#include "YinSet/SimplicialComplex.h"
#include "Core/Box.h"

template <int Order>
YinSet<2, Order> intersect(const YinSet<2, Order>& lhs,
                           const YinSet<2, Order>& rhs, Real tol);

template <int Dim, int Order>
struct CurveFactory;
template <int Order>
class YinSet<2, Order> : public SegmentedRealizableSpadjor<Order> {
 public:
  enum { Dim = 2 };
  using SRS = SegmentedRealizableSpadjor<Order>;
  using rVec = Vec<Real, 2>;
  using PointIndex = typename OrientedJordanCurve<Dim, Order>::PointIndex;
  using Vertex = typename Simplex<PointIndex>::Vertex;
  using YinSetPtr = std::shared_ptr<YinSet<Dim, Order>>;
  // using YinSetPtr = YinSet<Dim, Order>*;

  /// Initialize a Yin set from the stream.
  ///
  /// The Hasse diagram is calculated, but no re-segmentation is applied.
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
    int depth;  // even number for positive orientation, odd number of negative
                // orientation
    int parent;
    std::vector<int> children;
  };

  /// The Hasse diagram. The last node is the root of the forest.
  std::vector<Node> diagram;

 public:
  /// Get the boundary Jordan curves.
  const std::vector<OrientedJordanCurve<Dim, Order>>& getBoundaryCycles()
      const {
    return orientedJordanCurves;
  }

  /// Return if the Yin set is bounded, based on the Hasse diagram.
  bool isBounded() const { return (diagram.back().depth + 2) % 2 == 1; }

  /// Get the orientation of the k-th Jordan curve, based on the Hasse diagram.
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
  /// Calculate the complementary set.
  YinSet<2, Order> complement(Real tol) const;

  /// Calculate the intersection of two Yin sets.
  friend YinSet<2, Order> intersect<Order>(const YinSet<2, Order>& lhs,
                                           const YinSet<2, Order>& rhs,
                                           Real tol);

  // YinSet Factory
  friend struct CurveFactory<2, Order>;

  // kinks accessor
  const SimplicialComplex<Vertex>& getKinks() const { return kinks; }

  // initial kinks, remove old data. will refit every segmentedCurves.
  void resetAllKinks(std::vector<Vertex> vertices);

  // modify kinks, will refit related segmentedCurve.
  int insertKink(const Vertex& index);
  int eraseKink(const Vertex& index);

  // split orientedJordanCurves to spline without kinks.
  vector<Curve<2, Order>> getSmoothCurves(Real tol) const;

  /// Calculate CutCell.
  auto cutCell(const Box<Dim>& box, const Interval<Dim>& range,
               bool AddInner = false) const
      -> std::tuple<vector<vector<YinSetPtr>>,
                    vector<vector<vector<Curve<2, Order>>>>,
                    vector<vector<int>>>;

 protected:
  ///
  void buildHasse(Real tol);

  // refit segmentedCurve_i.
  void reFitCurve(size_t i);

  int bettiNumbers[2];
  SimplicialComplex<Vertex> kinks;
  using SRS::orientedJordanCurves;
};

template <int Order>
using YinSetPtr = YinSet<2, Order>::YinSetPtr;

template <int Order>
void dumpVecYinSet(const vector<vector<YinSetPtr<Order>>>& data,
                      std::ostream& os) {
  int num = 0;
  int N1 = data.size();
  int N2 = data[0].size();
  for (int i0 = 0; i0 < N1; i0++) {
    for (int i1 = 0; i1 < N2; i1++) {
      if (data[i0][i1]) {
        num++;
      }
    }
  }
  os.write((char*)&num, sizeof(num));
  for (int i0 = 0; i0 < N1; i0++) {
    for (int i1 = 0; i1 < N2; i1++) {
      if (data[i0][i1]) {
        data[i0][i1]->dump(os);
      }
    }
  }
}

#endif  // YINSET_H
