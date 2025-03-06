#ifndef SEGMENT_H
#define SEGMENT_H

#include "Core/Vec.h"

/// Define a line segment by two distinct endpoints
/**
 */
template <int Dim>
class Segment {
 public:
  typedef Vec<Real, Dim> rVec;

  rVec p[2];
  ///
  /**
     3 types of intersection for a pair of segments.
  */
  enum intsType { None = 0, One = 1, Overlap = 2 };

 public:
  ///
  /**
     The constructor accepting two endpoints.

     Segment(p1,p2) is equivalent to Segment(p2,p1).
   */
  Segment(const rVec &_p1, const rVec &_p2) {
    p[0] = _p1;
    p[1] = _p2;
  }

 public:
  ///
  /**
     Return the slope of a planar line segment. For 2-D only.

     May yield +inf -inf.
   */
  Real slope() const;

  ///
  /**
     Return the intersection point of the segment
     and the hyperplane x_{dim-1} = query.coord[dim-1]

     If the segment lies in that hyperplane, return query instead.
   */
  rVec substitute(const rVec &query, Real tol) const;

  ///
  /**
     Determine if two segments are parallel to each other.

     Two segments are considered parallel
     if the endpoints of their directional vectors (origins aligned)
     are about or less than the scale of tol.
   */
  bool parallel(const Segment<Dim> &that, Real tol) const;

  ///
  /**
     Return the signed area of the parallelogram formd by p[1]-p[0] and q-p[0].

     For 2-D only.
   */
  Real area(const rVec &q) const;

  ///
  /**
     Determine if the point q lies on the line segment.
   */
  bool contain(const rVec &q, Real tol) const;

  ///
  /**
   */
  friend std::ostream &operator<<(std::ostream &os, const Segment<Dim> &w) {
    os << w.p1 << "_" << w.p2;
    return os;
  }
};

//=================================================

template <>
inline Real Segment<2>::slope() const {
  rVec delta = p[0] - p[1];
  return delta[0] == 0
             ? std::numeric_limits<Real>::max() * (delta[1] > 0 ? 1 : -1)
             : (delta[1] / delta[0]);
}

template <int Dim>
inline Vec<Real, Dim> Segment<Dim>::substitute(const rVec &query,
                                               Real tol) const {
  Real dt = std::abs(p[0][Dim - 1] - p[1][Dim - 1]);
  // check if the query is valid
  if (std::abs(2 * query[Dim - 1] - p[0][Dim - 1] - p[1][Dim - 1]) > dt + tol)
    assert(!"Query out of range");
  if (dt <= tol)  // return query if the segment is contained in the hyperplane
    return query;

  Real ratio =
      (query[Dim - 1] - p[0][Dim - 1]) / (p[1][Dim - 1] - p[0][Dim - 1]);
  rVec r = p[0] + (p[1] - p[0]) * ratio;
  return r;
}

template <>
inline bool Segment<2>::parallel(const Segment<2> &that, Real tol) const {
  const rVec &p3 = that.p[0];
  const rVec &p4 = that.p[1];
  rVec A[2];
  A[0] = p[1] - p[0];
  A[1] = p3 - p4;
  return std::abs(cross(A[0], A[1])) / norm(A[0]) <= tol;
}

template <>
inline Real Segment<2>::area(const rVec &q) const {
  return cross(p[1] - p[0], q - p[0]);
}

// p1,p2 are endpoints of the first segment,
// while p3,p4 are that of the second.
template <int Dim>
inline typename Segment<Dim>::intsType solveForOverlie(
    Vec<Real, Dim> &p1, Vec<Real, Dim> &p2, Vec<Real, Dim> &p3,
    Vec<Real, Dim> &p4, Vec<Real, Dim> &pOut1, Vec<Real, Dim> &pOut2, Real tol,
    int majorDim) {
  if (p1[majorDim] > p2[majorDim]) std::swap(p1, p2);
  if (p3[majorDim] > p4[majorDim]) std::swap(p3, p4);
  pOut1 = (p1[majorDim] < p3[majorDim]) ? (p3) : (p1);
  pOut2 = (p2[majorDim] < p4[majorDim]) ? (p2) : (p4);
  Real r = pOut2[majorDim] - pOut1[majorDim];
  if (r < -tol) {
    return Segment<Dim>::intsType::None;
  } else if (r > tol) {
    return Segment<Dim>::intsType::Overlap;
  }
  pOut1 = (pOut1 + pOut2) * 0.5;
  return Segment<Dim>::intsType::One;
}

inline Segment<2>::intsType intersect(const Segment<2> &thisSeg,
                                      const Segment<2> &thatSeg,
                                      Vec<Real, 2> &pOut1, Vec<Real, 2> &pOut2,
                                      Real tol) {
  using rVec = Segment<2>::rVec;
  rVec p1 = thisSeg.p[0];
  rVec p2 = thisSeg.p[1];
  rVec p3 = thatSeg.p[0];
  rVec p4 = thatSeg.p[1];

  rVec A[2], b;
  A[0] = p2 - p1;
  A[1] = p3 - p4;
  b = p3 - p1;
  Real det = cross(A[0], A[1]);
  Real sc = norm(A[0]);

  if (std::abs(det / sc) < tol) {  // parallel segments
    Real r = cross(A[0], b) / sc;
    if (std::abs(r) > tol) return Segment<2>::intsType::None;
    int majorDim = (std::abs(A[0][0]) > std::abs(A[0][1])) ? (0) : (1);
    return solveForOverlie(p1, p2, p3, p4, pOut1, pOut2, tol, majorDim);
  }

  // solve for intersections by Cramer's rule
  Real x[2];
  x[0] = cross(b, A[1]) / det;
  x[1] = cross(A[0], b) / det;

  if (x[0] > -tol / sc && x[0] < 1 + tol / sc && x[1] > -tol / sc &&
      x[1] < 1 + tol / sc) {
    pOut1 = p1 + (p2 - p1) * x[0];
    return Segment<2>::intsType::One;
  }
  return Segment<2>::intsType::None;
}

template <int Dim>
inline bool isInInterval(const Vec<Real, Dim> &q, const Vec<Real, Dim> &p1,
                         const Vec<Real, Dim> &p2, Real tol, int majorDim) {
  return static_cast<bool>(q[majorDim] >= std::min(p1[majorDim], p2[majorDim]) - tol &&
      q[majorDim] <= std::max(p1[majorDim], p2[majorDim]) + tol);
}

template <int Dim>
inline bool Segment<Dim>::contain(const rVec &q, Real tol) const {
  if (!parallel(Segment<Dim>(q, p[0]), tol)) return false;
  auto delta = p[1] - p[0];
  int majorDim = (std::abs(delta[0]) > std::abs(delta[1])) ? (0) : (1);
  return isInInterval(q, p[0], p[1], tol, majorDim);
}

#endif  // SEGMENT_H
