#ifndef SEGMENTSINTERSECTOR_H
#define SEGMENTSINTERSECTOR_H

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "Core/VecCompare.h"
#include "Segment.h"

/// Compute the intersections of a multiset of line segments.
/**
All the intersection points along with the incident line segments will be
reported. Particularly, for an overlie, endpoints of the overlie will be
reported.
 */

class SegmentsIntersector {
 public:
  typedef Vec<Real, 2> rVec;
  typedef Segment<2> Seg;
  typedef typename std::vector<Seg>::const_iterator cpSeg;
  typedef std::vector<cpSeg> Bundle;

  typedef std::map<rVec, std::vector<int>, VecCompare<Real, 2>> ResultType;

 protected:
  ///
  /**
     A partial order defined on bundles
     is required by the line sweeping algorithm.

     The x-coordinate of the intersection of the sweep line and the bundle
     is the major keyword.
     The slope is the minor keyword.
   */
  struct BundleCompare {
    Real tol;
    const rVec *volatile evt;  // the event point, volatile ?

    BundleCompare() : tol(0), evt(nullptr) {}
    BundleCompare(Real _tol, const rVec *_q) : tol(_tol), evt(_q) {}

    bool operator()(const Bundle *b1, const Bundle *b2) const {
      assert(!b1->empty());
      assert(!b2->empty());
      // find representative from b1, b2
      cpSeg seg1 = b1->front();
      cpSeg seg2 = b2->front();
      Real x1 = seg1->substitute(*evt, tol)[0];
      Real x2 = seg2->substitute(*evt, tol)[0];
      if (std::abs(x1 - x2) < tol) {  // use slope
        if (seg1->parallel(*seg2, tol)) return false;
        Real k1 = seg1->slope();
        Real k2 = seg2->slope();
        if (k1 == 0 || k2 == 0) return (k1 != 0);
        if (k1 * k2 < 0) return (k1 > 0);
        return (k1 < k2);
      }
      return x1 < x2;
    }
  };

  typedef std::map<rVec, std::map<cpSeg, int>, VecCompare<Real, 2>> QueueType;
  typedef std::set<Bundle *, BundleCompare> StatusType;

 public:
  ///
  /**
     The constructor accepts _tol as the uncertainty quantification.
   */
  SegmentsIntersector(Real _tol)
      : tol(_tol), pt_cmp(_tol), bdl_cmp(_tol, &eventPoint) {}

  ///
  /**
     The copy constructor is disabled.
   */
  SegmentsIntersector(const SegmentsIntersector &) = delete;

  ///
  /**
     This function does the work of computing the intersections.

     The uncertainty amount used here has been passed in through the
     constructor.

     Return a std::map<> object (see ResultType) whose keys are the intersection
     points, and whose values are the incident segments (represented by their
     indices in lineSet).
   */
  ResultType operator()(const std::vector<Seg> &lineSet) {
    event.clear();
    event = QueueType(pt_cmp);
    status.clear();
    status = StatusType(bdl_cmp);

    lineSetBegin = lineSet.cbegin();
    backref.resize(lineSet.size());
    return findIntersections(lineSet);
  }

 protected:
  ///
  /**
     The main function for finding the intersections
     for a multiset of line segments.

     Refer to findIntersections() in the document for details.
   */
  ResultType findIntersections(const std::vector<Seg> &lineSet);

  ///
  /**
     Insert an event point with one incident segment
     into the queue.
   */
  void eventEnqueue(const rVec &p, cpSeg cit, int type);

  ///
  /**
     For each pair of segments from the two bundles respectively,
     test for intersections.

     Refer to findNewEvents() in the document for details.
   */
  void findNewEvents(const Bundle *s1, const Bundle *s2);

  ///
  /**
     Test a pair of segments for intersections.
     This is used by the bundle version of findNewEvents().
   */
  void findNewEvents(cpSeg seg1, cpSeg seg2);

  ///
  /**
     Recover the true C(p) under special case.

     Refer to the document for detailed dicussion on the algorithm.
   */
  void collectCp(cpSeg testor);

  ///
  /**
     Add a segment to the status structure.

     If there is an existing bundle to which the segment belongs to,
     the segment will be added to that bundle
     and tested against all the other segments in that bundle.
     Otherwise a new bundle is created, containing only the segment.

     Refer to addToStatus() in the document for details.
   */
  StatusType::iterator addToStatus(cpSeg cit);

  ///
  /**
     Delete a segment from the status structure.
   */
  StatusType::iterator delFromStatus(cpSeg cit);

 protected:
  enum { INTERIOR = 0, UPPER = -1, LOWER = 1 };
  std::vector<Seg>::const_iterator lineSetBegin;
  Real tol;

  rVec eventPoint;
  std::vector<StatusType::iterator> backref;
  QueueType event;
  StatusType status;
  VecCompare<Real, 2> pt_cmp;
  BundleCompare bdl_cmp;
};

#endif  // SEGMENTSINTERSECTOR_H
