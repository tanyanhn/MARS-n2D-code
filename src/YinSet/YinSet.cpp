#include "YinSet.h"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <queue>
#include <sstream>
#include "Core/Tensor.h"
#include "PastingMap.h"
#include "PointsLocater.h"

// re-declarations; their definitions are in SegmentedRealizableSpadjor.cpp.
template <int Order> std::vector<Segment<2>> collapseToSeg(const Curve<2, Order> &polygon);
template <int Order> bool isBounded(const Curve<2, Order> &polygon, Real tol);

//============================================================

template <int Order>
YinSet<2, Order>::YinSet(std::istream &is, Real tol)
    : SRS(is, 0.0)
{
  buildHasse(tol);
}

// template<> YinSet<2,2>::YinSet(std::istream &is, Real tol);
// template<> YinSet<2,4>::YinSet(std::istream &is, Real tol);

template <int Order>
YinSet<2, Order>::YinSet(const SRS& segmentedSpadjor, Real tol) {
  this->segmentedCurves = segmentedSpadjor.segmentedCurves;
  buildHasse(tol);
}

// template YinSet<2,2>::YinSet(const SRS &segmentedSpadjor, Real tol);
// template YinSet<2, 4>::YinSet(const SRS& segmentedSpadjor, Real tol);

template <int Order>
void YinSet<2, Order>::buildHasse(Real tol)
{
  // step 1 : construct the inclusion matrix
  const int numCurves = segmentedCurves.size();
  if (numCurves == 0) {
    bettiNumbers[0] = 0;
    bettiNumbers[1] = 0;
    return;
  }
  Tensor<int, 2> mat(Vec<int, 2>{numCurves, numCurves});
  mat = 0;
  std::vector<Vec<Real, 2>> somePoints(numCurves);
  std::vector<Interval<2>> boxes(numCurves);
  std::vector<bool> boundedness(numCurves);
  std::vector<std::vector<int>> candidates(numCurves);
  for(int i=0; i<numCurves; ++i) {
    const auto &polys = segmentedCurves[i].getPolys();
    const auto &knots = segmentedCurves[i].getKnots();
    somePoints[i] = polys[0]((knots[0] + knots[1]) / 2);
    boxes[i] = boundingBox(segmentedCurves[i]);
    boundedness[i] = ::isBounded(segmentedCurves[i], tol);
  }
  for(int i=0; i<numCurves-1; ++i) {
    for(int j=i+1; j<numCurves; ++j) {
      if(boxes[i].contain(boxes[j], tol))
        candidates[i].push_back(j);
      if(boxes[j].contain(boxes[i], tol))
        candidates[j].push_back(i);
    }
  }
  //
  PointsLocater locater(tol);
  for(int i=0; i<numCurves; ++i) {
    std::vector<Vec<Real,2>> queries;
    for(int k : candidates[i])
      queries.push_back(somePoints[k]);
    auto loc = locater(collapseToSeg(segmentedCurves[i]), queries, boundedness[i]);
    for(std::size_t j = 0; j < candidates[i].size(); ++j) {
      assert(loc[j] != 0);
      if(boundedness[i] == (loc[j] == 1)) {
        mat(i, candidates[i][j]) = 1;
        mat(candidates[i][j], i) = -1;
      }
    }
  }
  // step 2 : obtain the Hasse diagram from the inclusion matrix
  diagram.resize(numCurves+1);
  diagram[numCurves].depth = -2;
  diagram[numCurves].parent = -1;
  std::vector<int> numAnc(numCurves); // number of ancestors
  std::vector<int> parent(numCurves);
  for(int i = 0; i < numCurves; ++i) {
    numAnc[i] = std::count(&mat(0, i), &mat(0, i+1), 1);
    if(numAnc[i] == 0) {
      parent[i] = numCurves;
      diagram[numCurves].depth = (boundedness[i]) ? (-1) : (0);
    }
  }
  assert(diagram[numCurves].depth != -2);
  // topological sort
  int numRemain = numCurves;
  int numPositive = 0;
  while(numRemain--) {
    int j = std::find(numAnc.cbegin(), numAnc.cend(), 0) - numAnc.cbegin();
    numAnc[j] = -1; // so that j will never be selected again
    diagram[j].parent = parent[j];
    diagram[j].depth = diagram[parent[j]].depth + 1;
    if(diagram[j].depth % 2 == 0)
      ++numPositive;
    diagram[parent[j]].children.push_back(j);
    // handle the children
    for(int k = 0; k < numCurves; ++k) {
      if(mat(j, k) == 1) {
        parent[k] = j;
        --numAnc[k];
      }
    }
  }
  // step 3 : calculate the Betti nubmers.
  if(diagram[numCurves].depth == -1) {
    // this Yin set is bounded.
    bettiNumbers[0] = numPositive;
    bettiNumbers[1] = numCurves - numPositive;
  } else {
    // this Yin set is unbounded.
    bettiNumbers[0] = numPositive + 1;
    bettiNumbers[1] = numCurves - numPositive;
  }
}

//============================================================

template <int Order>
YinSet<2, Order> YinSet<2, Order>::complement(Real tol) const {
  return YinSet<2, Order>(SegmentedRealizableSpadjor<Order>::complement(tol),
                          tol);
}

template <>
YinSet<2, 2> intersect(const YinSet<2, 2>& lhs,
                       const YinSet<2, 2>& rhs,
                       Real tol) {
  auto segmented = meet(lhs, rhs, tol);\
  return YinSet<2, 2>(segmented, tol);
}

//============================================================

bool equal(const Curve<2,2> &lhs, const Curve<2,2> &rhs, Real tol)
{
  std::size_t N;
  if((N = lhs.getPolys().size()) != rhs.getPolys().size())
    return false;
  VecCompare<Real,2> vcmp(tol);
  auto findTopmostIdx = [&vcmp](const Curve<2,2> &G) {
    Vec<Real,2> r(-std::numeric_limits<Real>::max());
    std::size_t a;
    const auto &polys = G.getPolys();
    for(std::size_t i=0; i<polys.size(); ++i)
      if(vcmp.compare(polys[i][0], r) == -1) {
        r = polys[i][0];
        a = i;
      }
    return a;
  };
  std::size_t k1 = findTopmostIdx(lhs);
  std::size_t k2 = findTopmostIdx(rhs);
  for(std::size_t i=0; i<N; ++i) {
    if(vcmp(lhs.getPolys()[k1][0], rhs.getPolys()[k2][0]) != 0)
      return false;
    k1 = (k1+1) % N;
    k2 = (k2+1) % N;
  }
  return true;
}

template <>
bool YinSet<2, 2>::equal(const YinSet<2, 2> &rhs, Real tol) const
{
  VecCompare<Real,2> vcmp(tol);
  using Item = std::pair<int, bool>;
  std::map<Vec<Real,2>, std::vector<Item>, VecCompare<Real,2>> lookup(vcmp);

  auto findTopmost = [&vcmp](const Curve<2,2> &G) {
    Vec<Real,2> r(-std::numeric_limits<Real>::max());
    const auto &polys = G.getPolys();
    for(std::size_t i = 0; i < polys.size(); ++i)
      if(vcmp.compare(polys[i][0], r) == -1)
        r = polys[i][0];
    return r;
  };
  for(std::size_t i = 0; i < segmentedCurves.size(); ++i) {
    Vec<Real,2> topmost = findTopmost(segmentedCurves[i]);
    auto ret = lookup.insert(std::make_pair(topmost, std::vector<Item>()));
    ret.first->second.push_back(std::make_pair(i, false));
  }
  for(const auto &j : rhs.segmentedCurves) {
    Vec<Real,2> topmost = findTopmost(j);
    auto iter = lookup.find(topmost);
    if(iter == lookup.end())
      return false;
    std::vector<Item> &items = iter->second;
    auto k = items.begin();
    for(; k!=items.end(); ++k) {
      if(!k->second && ::equal(j, segmentedCurves[k->first], tol)) {
        k->second = true;
        break;
      }
    }
    if(k == items.end())
      return false;
  }
  return true;
}

//============================================================

template <int Order>
std::string YinSet<2, Order>::getHasseString() const
{
  if(segmentedCurves.empty()) return "YinSet empty.";
  const int w = 8;
  std::ostringstream oss;
  oss << std::left << std::setw(w) << " ";
  oss << std::left << std::setw(w) << "Parent";
  oss << std::left << std::setw(w) << "Depth";
  oss << std::left << std::setw(w) << "Orient";
  oss << std::left << "Children" << "\n";

  for(std::size_t i=0; i<diagram.size(); ++i) {
    oss << std::left << std::setw(w) << i;
    oss << std::left << std::setw(w) << diagram[i].parent;
    oss << std::left << std::setw(w) << diagram[i].depth;
    if(i < segmentedCurves.size()) {
      oss << std::left << std::setw(w) << getOrientation(i);
    } else {
      oss << std::left << std::setw(w) << " ";
    }
    oss << "{";
    const std::vector<int> &children = diagram[i].children;
    for(int j : children)
      oss << j << ",";
    oss << "}" << "\n";
  }
  return oss.str();
}

template <int Order>
void YinSet<2, Order>::dump(std::ostream &os) const
{
  const int N = segmentedCurves.size();
  os.write((char*)&N, sizeof(int));
  // save the boundaries according to the order of BFS
  const auto &rootOfForest = diagram.back();
  std::queue<int> Q;
  for(int i : rootOfForest.children)
    Q.push(i);
  while(!Q.empty()) {
    int k = Q.front();
    Q.pop();
    for(int i : diagram[k].children)
      Q.push(i);
    segmentedCurves[k].dump(os);
  }
}

template <int Order>
void YinSet<2, Order>::setKinks(
    std::vector<std::pair<unsigned int, unsigned int>> vertices) {
  sort(vertices.begin(), vertices.end());
  mPoint2Vertex.clear();
  mVertex2Point.clear();
  kinks = SimplicialComplex();
  unsigned int vertex = 0;
  auto start = vertices.begin(), end = start;
  auto numCurves = segmentedCurves.size();
  for (size_t i = 0; i < numCurves; ++i) {
    start = std::lower_bound(vertices.begin(), vertices.end(),
                             std::make_pair<unsigned int, unsigned int>(i, 0));
    end =
        std::lower_bound(vertices.begin(), vertices.end(),
                         std::make_pair<unsigned int, unsigned int>(i + 1, 0));
    for (auto Iter = start; Iter != end; ++Iter) {
      mPoint2Vertex[*Iter] = vertex;
      mVertex2Point[vertex] = *Iter;
      kinks.insert(Simplex{std::initializer_list<unsigned int>{vertex}});
      ++vertex;
    }
    reFitCurve(i);
  }
}

template <int Order>
int YinSet<2, Order>::vertex2Point(
    unsigned int vertex,
    std::pair<unsigned int, unsigned int>& index) const {
  auto iter = mVertex2Point.find(vertex);
  int ret = 1;
  if (iter == mVertex2Point.end()) {
    ret = 0;
  } else {
    index = iter->second;
  }
  return ret;
}

template <int Order>
int YinSet<2, Order>::vertex2Point(unsigned int vertex, rVec& point) const {
  std::pair<unsigned int, unsigned int> index;
  int ret = vertex2Point(vertex, index);
  if (ret == 1)
    point = segmentedCurves[index.first](
        segmentedCurves[index.first].getKnots()[index.second]);
  return ret;
}

template <int Order>
int YinSet<2, Order>::point2Vertex(
    const std::pair<unsigned int, unsigned int>& index,
    unsigned int& vertex) const {
  auto iter = mPoint2Vertex.find(index);
  int ret = 1;
  if (iter == mPoint2Vertex.end()) {
    ret = 0;
  } else {
    vertex = iter->second;
  }
  return ret;
}

template <int Order>
int YinSet<2, Order>::insertKinks(
    const std::pair<unsigned int, unsigned int>& index) {
  if (mPoint2Vertex.find(index) != mPoint2Vertex.end())
    return -1;
  auto i = index.first;
  unsigned int vertex;
  if (mVertex2Point.empty())
    vertex = 0;
  else
    vertex = mVertex2Point.rbegin()->first + 1;
  mPoint2Vertex[index] = vertex;
  mVertex2Point[vertex] = index;
  kinks.insert(Simplex{std::initializer_list<unsigned int>{vertex}});
  reFitCurve(i);
  return vertex;
}

template <int Order>
int YinSet<2, Order>::eraseKinks(unsigned int vertex) {
  if (mVertex2Point.find(vertex) == mVertex2Point.end())
    return -1;
  auto& index = mVertex2Point[vertex];
  auto i = index.first;
  mPoint2Vertex.erase(index);
  mVertex2Point.erase(vertex);
  kinks.erase(Simplex{std::initializer_list<unsigned int>{vertex}});
  reFitCurve(i);
  return vertex;
}

template <int Order>
void YinSet<2, Order>::reFitCurve(unsigned int i) {
  auto start = mPoint2Vertex.lower_bound(std::make_pair(i, 0)),
       end = mPoint2Vertex.lower_bound(std::make_pair(i + 1, 0));
  SimplicialComplex sims;
  vector<Vec<Real, 2>> points;
  while (start != end) {
    sims.insert(
        Simplex{std::initializer_list<unsigned int>{start->first.second}});
    start++;
  }
  for (auto t : segmentedCurves[i].getKnots())
    points.push_back(segmentedCurves[i](t));
  segmentedCurves[i].define(points, sims);
}

//============================================================
template class YinSet<2, 2>;
template class YinSet<2, 4>;
