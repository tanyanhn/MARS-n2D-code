#ifndef POINTSLOCATER_H
#define POINTSLOCATER_H

#include <set>
#include <vector>

#include "Core/VecCompare.h"
#include "Segment.h"
#include "YinSet/OrientedJordanCurve.h"

/// Determine the relative locations of a set of points to a spadjor forest / an
/// oriented polygon
/**
   Refer to the Algorithm locate().
 */
class PointsLocater {
 public:
  using rVec = Vec<Real, 2>;
  template <class T>
  using vector = std::vector<T>;

  ///
  /**
     The constructor accepts the absolute tolerance.
   */
  PointsLocater(Real _tol) : tol(_tol), pt_cmp(_tol) {
    seg_cmp = SegCompare(_tol, &ep);
    status = StatusType(seg_cmp);
  }

  ///
  /**
   */
  PointsLocater(const PointsLocater &) = delete;

  vector<int> operator()(const vector<Segment<2>> &ys,
                         const vector<rVec> &queries, bool bounded) {
    return compute(ys, queries, (bounded) ? (-1) : (1));
  }
  template <int Order>
  vector<int> operator()(
      const vector<OrientedJordanCurve<DIM, Order>> &ys,
      const vector<rVec> &queries,
      bool bounded);

 protected:
  vector<int> compute(const vector<Segment<2>> &segs, const vector<rVec> &pts,
                      int defVal);

  using cpSeg = typename vector<Segment<2>>::const_iterator;

  struct SegCompare {
    Real tol;
    const rVec *volatile query;

    SegCompare() : tol(0), query(0) {}
    SegCompare(Real _tol, const rVec *_q) : tol(_tol), query(_q) {}

    bool operator()(cpSeg s1, cpSeg s2) const {
      Real x1 = s1->substitute(*query, tol)[0];
      Real x2 = s2->substitute(*query, tol)[0];
      if (std::abs(x1 - x2) < tol) {  // use slope
        if (s1->parallel(*s2, tol)) return false;
        Real k1 = s1->slope();
        Real k2 = s2->slope();
        if (k1 == 0 || k2 == 0) return (k1 != 0);
        if (k1 * k2 < 0) return (k1 > 0);
        return (k1 < k2);
      }
      return x1 < x2;
    }
  };

  struct Event {
    rVec p;
    int type;  // -1=upper, -3=query, 1=lower
    int outIdx;
    cpSeg cit;
    Event(const rVec &_p, int _type, int _outIdx, cpSeg _cit)
        : p(_p), type(_type), outIdx(_outIdx), cit(_cit) {}
  };

  using StatusType = std::set<cpSeg, SegCompare>;

  Real tol;
  rVec ep;
  std::vector<StatusType::const_iterator> backref;
  std::vector<Event> event;
  StatusType status;

  VecCompare<Real, 2> pt_cmp;
  SegCompare seg_cmp;
};

#endif  // POINTSLOCATER_H
