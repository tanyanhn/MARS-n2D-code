#include "SegmentedRealizableSpadjor.h"
#include <limits>
#include "PointsLocater.h"
#include "SegmentsIntersector.h"

template <class T>
bool checkIsBounded(const T& iterator, Real tol) {
  VecCompare<Real, 2> vcmp(tol);
  Vec<Real, 2> topmost(-std::numeric_limits<Real>::max());
  std::vector<Segment<2>> edges;
  // collect the segments incident to the topmost vertex
  iterator([&](const Segment<2>& seg) {
    for (int i = 0; i < 2; ++i) {
      auto r = vcmp.compare(seg.p[i], topmost);
      if (r == -1) {
        topmost = seg.p[i];
        edges.clear();
        edges.push_back(seg);
      } else if (r == 0) {
        edges.push_back(seg);
      }
    }
  });
  // find the leftmost segment among the ones incident to the topmost vertex
  std::sort(edges.begin(), edges.end(),
            [](const Segment<2>& lhs, const Segment<2>& rhs) {
              Real k1 = lhs.slope();
              Real k2 = rhs.slope();
              if (k1 == 0)
                return false;
              if (k2 == 0)
                return true;
              if (k1 * k2 < 0)
                return k1 > 0;
              return k1 < k2;
            });
  return vcmp.compare(topmost, edges[0].p[0]) == 0;
}

template <>
bool SegmentedRealizableSpadjor<2>::isBounded(Real tol) const {
  auto iterator = [&](const auto& callback) {
    for (const auto& gamma : segmentedCurves) {
      const auto& polys = gamma.getPolys();
      int numPolys = polys.size();
      for (int i = 0; i < numPolys - 1; ++i)
        callback(Segment<2>(polys[i][0], polys[i + 1][0]));
      callback(Segment<2>(polys[numPolys - 1][0], polys[0][0]));
    }
  };
  return checkIsBounded(iterator, tol);
}

template <int Order>
bool isBounded(const Curve<2, Order>& polygon, Real tol) {
  auto iterator = [&](const auto& callback) {
    const auto& polys = polygon.getPolys();
    int numPolys = polys.size();
    for (int i = 0; i < numPolys - 1; ++i)
      callback(Segment<2>(polys[i][0], polys[i + 1][0]));
    callback(Segment<2>(polys[numPolys - 1][0], polys[0][0]));
  };
  return checkIsBounded(iterator, tol);
}

template bool isBounded(const Curve<2, 2>&, Real);
template bool isBounded(const Curve<2, 4>&, Real);

//============================================================

template <int Order>
std::vector<Segment<2>> collapseToSeg(const Curve<2, Order>& polygon) {
  std::vector<Segment<2>> segs;
  const auto& polys = polygon.getPolys();
  auto np = polys.size();
  for (std::size_t i = 0; i < np - 1; ++i)
    segs.emplace_back(polys[i][0], polys[i + 1][0]);
  segs.emplace_back(polys[np - 1][0], polygon.endpoint());
  return segs;
}

template std::vector<Segment<2>> collapseToSeg(const Curve<2, 2>&);
template std::vector<Segment<2>> collapseToSeg(const Curve<2, 4>&);

std::vector<Segment<2>> collapseToSeg(const std::vector<Curve<2, 2>>& bdries,
                                      std::vector<int>& cIdx,
                                      std::vector<int>& pIdx) {
  std::vector<Segment<2>> segs;
  for (std::size_t c = 0; c < bdries.size(); ++c) {
    const auto& gamma = bdries[c];
    const auto& polys = gamma.getPolys();
    auto numPolys = polys.size();
    for (std::size_t i = 0; i < numPolys - 1; ++i) {
      segs.emplace_back(polys[i][0], polys[i + 1][0]);
      cIdx.push_back(c);
      pIdx.push_back(i);
    }
    segs.emplace_back(polys[numPolys - 1][0], gamma.endpoint());
    cIdx.push_back(c);
    pIdx.push_back(numPolys - 1);
  }
  return segs;
}

template <>
void SegmentedRealizableSpadjor<2>::applySegmentation(
    const vector<Curve<Dim, 2>>& aSpadjor,
    Real tolForSegmentation) {
  std::vector<int> idxOfCurve;
  std::vector<int> idxOfPiece;
  std::vector<Segment<2>> segs =
      collapseToSeg(aSpadjor, idxOfCurve, idxOfPiece);
  SegmentsIntersector intersector(tolForSegmentation);
  auto itsInfo = intersector(segs);
  // Find the curve parameters at the improper intersections
  std::vector<std::vector<Real>> brks(
      aSpadjor.size());  // [which curve][which piece]
  for (const auto& info : itsInfo) {
    const auto& incident = info.second;
    if (incident.size() > 2) {
      for (int j : incident) {
        const auto& knots = aSpadjor[idxOfCurve[j]].getKnots();
        const auto& polys = aSpadjor[idxOfCurve[j]].getPolys();
        auto lastknot = polys[idxOfPiece[j]][0];
        Real t = knots[idxOfPiece[j]] + norm(info.first - lastknot, 2);
        brks[idxOfCurve[j]].push_back(t);
      }
    } else {
      assert(incident.size() == 2);
      if (incident[0] != incident[1])
        assert(
            !"Invalid input of Yin set : may contain proper intersections. ");
    }
  }
  // split the curves
  for (std::size_t k = 0; k < aSpadjor.size(); ++k) {
    std::sort(brks[k].begin(), brks[k].end());
    aSpadjor[k].split(brks[k], segmentedCurves, tolForSegmentation);
  }
}

template <>
SegmentedRealizableSpadjor<2> meet(const SegmentedRealizableSpadjor<2>& lhs,
                                   const SegmentedRealizableSpadjor<2>& rhs,
                                   Real tol) {
  const int Dim = 2;
  const int Order = 2;
  using std::pair;
  using std::vector;
  const SegmentedRealizableSpadjor<Order>* operands[] = {&lhs, &rhs};
  const std::vector<Curve<Dim, Order>>* bdriesOfOperands[] = {
      &(lhs.segmentedCurves), &(rhs.segmentedCurves)};
  SegmentedRealizableSpadjor<Order> result;
  // group all the line segments
  vector<Segment<Dim>> segs[2];
  vector<int> ci;  // index of curve
  vector<int> pi;  // index of piece
  int numSeg[2];
  for (int m = 0; m < 2; ++m) {
    segs[m] = collapseToSeg(*bdriesOfOperands[m], ci, pi);
    numSeg[m] = segs[m].size();
  }

  // call the intersector
  vector<pair<Vec<Real, Dim>, vector<int>>> itsInfo;
  {
    vector<Segment<2>> S;
    S.insert(S.cend(), segs[0].cbegin(), segs[0].cend());
    S.insert(S.cend(), segs[1].cbegin(), segs[1].cend());
    SegmentsIntersector its(tol);
    auto rawInfo = its(S);
    for (auto& p : rawInfo) {
      const auto& incident = p.second;
      if (incident.size() > 2 ||
          (incident[0] < numSeg[0]) != (incident[1] < numSeg[0]))
        itsInfo.push_back(std::make_pair(p.first, std::move(p.second)));
    }
  }

  // dispatch the intersections to each curve
  vector<vector<Real>> breaks[2];  // [which SF][which curve][which break]
  for (int m = 0; m < 2; m++)
    breaks[m].resize(bdriesOfOperands[m]->size());

  for (const auto& i : itsInfo) {
    for (int j : i.second) {
      // for every Jordan curve in the realizable spadjor, find out all the
      // breaks
      int m = (j < numSeg[0]) ? (0) : (1);
      // find the parameter
      const auto& knots = (*bdriesOfOperands[m])[ci[j]].getKnots();
      const auto& polys = (*bdriesOfOperands[m])[ci[j]].getPolys();
      auto lastknot = polys[pi[j]][0];
      Real t = knots[pi[j]] + norm(i.first - lastknot);
      breaks[m][ci[j]].push_back(t);
    }
  }

  // extract the beta_i
  vector<Curve<Dim, Order>> beta[2];
  for (int m = 0; m < 2; ++m) {
    for (std::size_t c = 0; c < bdriesOfOperands[m]->size(); ++c) {
      // sort the breaks
      std::sort(breaks[m][c].begin(), breaks[m][c].end());
      // split() deals with tiny pieces
      (*bdriesOfOperands[m])[c].split(breaks[m][c], beta[m], tol);
    }
  }

  vector<Vec<Real, Dim>> midBeta[2];
  vector<int> locs[2];
  PointsLocater locater(tol);
  // determine which beta_i is in the interior of the other spadjor forest
  for (int m = 0; m < 2; ++m) {
    midBeta[m].resize(beta[m].size());
    for (std::size_t i = 0; i < beta[m].size(); ++i)
      midBeta[m][i] = beta[m][i].midpoint();
    locs[m] = locater(segs[1 - m], midBeta[m], operands[1 - m]->isBounded(tol));
    for (std::size_t i = 0; i < beta[m].size(); ++i)
      if (locs[m][i] == 1)
        result.segmentedCurves.push_back(std::move(beta[m][i]));
  }

  // deal with the overlapping edges
  VecCompare<Real, Dim> vcmp(tol);
  std::map<Vec<Real, Dim>, int, VecCompare<Real, Dim>> mid2idx(vcmp);
  for (std::size_t i = 0; i < midBeta[1].size(); ++i) {
    if (locs[1][i] == 0)
      mid2idx.insert(std::make_pair(midBeta[1][i], i));
  }
  for (std::size_t i = 0; i < beta[0].size(); ++i) {
    if (locs[0][i] == 0) {
      auto q = midBeta[0][i];
      assert(mid2idx.find(q) != mid2idx.end());
      int k = mid2idx[q];
      // is of same orientation ?
      if (norm(beta[0][i].startpoint() - beta[1][k].startpoint()) < tol)
        result.segmentedCurves.push_back(std::move(beta[0][i]));
    }
  }
  return result;
}

//============================================================

template <>
SegmentedRealizableSpadjor<2>::SegmentedRealizableSpadjor(
    const vector<Curve<Dim, 2>>& aSpadjor,
    Real tolForSegmentation) {
  if (tolForSegmentation == 0.0) {
    segmentedCurves = aSpadjor;
  } else {
    applySegmentation(aSpadjor, tolForSegmentation);
  }
}

template <>
SegmentedRealizableSpadjor<4>::SegmentedRealizableSpadjor(
    const vector<Curve<Dim, 4>>& aSpadjor,
    Real tolForSegmentation) {
  // *******************************************************
  // The implementation of the Boolean algebra for Yin sets
  // represented by cubic splines shoud go here.
  // assert(tolForSegmentation == 0.0);
  // *******************************************************
  segmentedCurves = aSpadjor;
}

template <>
SegmentedRealizableSpadjor<2>::SegmentedRealizableSpadjor(
    std::istream& is,
    Real tolForSegmentation) {
  vector<Curve<Dim, 2>> content;
  int N;
  is.read((char*)&N, sizeof(int));
  for (int i = 0; i < N; ++i)
    content.push_back(Curve<Dim, 2>::load(is));
  if (tolForSegmentation == 0.0) {
    segmentedCurves = std::move(content);
  } else {
    applySegmentation(content, tolForSegmentation);
  }
}

template <>
SegmentedRealizableSpadjor<4>::SegmentedRealizableSpadjor(
    std::istream& is,
    Real tolForSegmentation) {
  vector<Curve<Dim, 4>> content;
  int N;
  is.read((char*)&N, sizeof(int));
  for (int i = 0; i < N; ++i)
    content.push_back(Curve<Dim, 4>::load(is));
  // *******************************************************
  // The implementation of the Boolean algebra for Yin sets
  // represented by cubic splines shoud go here.
  assert(tolForSegmentation == 0.0);
  segmentedCurves = std::move(content);
  // *******************************************************
}

template <int Order>
auto SegmentedRealizableSpadjor<Order>::translate(const rVec& delta) const
    -> SRS {
  SRS result;
  for (const auto& j : segmentedCurves)
    result.segmentedCurves.push_back(j.offset(delta));
  return result;
}

template <int Order>
auto SegmentedRealizableSpadjor<Order>::complement() const -> SRS {
  SRS result;
  for (const auto& j : segmentedCurves)
    result.segmentedCurves.push_back(j.reverse());
  return result;
}

//============================================================
template class SegmentedRealizableSpadjor<2>;
template class SegmentedRealizableSpadjor<4>;