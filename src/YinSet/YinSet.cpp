#include "YinSet.h"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <queue>
#include <sstream>

#include "Core/Tensor.h"
#include "PointsLocater.h"
#include "YinSet/CutCellHelper.h"
// #include "YinSet/PastingMap.h"

// re-declarations; their definitions are in SegmentedRealizableSpadjor.cpp.
template <int Order>
std::vector<Segment<2>> collapseToSeg(const Curve<2, Order> &polygon);
template <int Order>
bool isBounded(const Curve<2, Order> &polygon, Real tol);

//============================================================

template <int Order>
YinSet<2, Order>::YinSet(std::istream &is, Real tol) : SRS(is, 0.0) {
  buildHasse(tol);
}

// template<> YinSet<2,2>::YinSet(std::istream &is, Real tol);
// template<> YinSet<2,4>::YinSet(std::istream &is, Real tol);

template <int Order>
YinSet<2, Order>::YinSet(const SRS &segmentedSpadjor, Real tol) {
  this->orientedJordanCurves = segmentedSpadjor.orientedJordanCurves;
  buildHasse(tol);
}

// template YinSet<2,2>::YinSet(const SRS &segmentedSpadjor, Real tol);
// template YinSet<2, 4>::YinSet(const SRS& segmentedSpadjor, Real tol);

template <int Order>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
void YinSet<2, Order>::buildHasse(Real tol) {
  // step 1 : construct the inclusion matrix
  const int numCurves = orientedJordanCurves.size();
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
  for (int i = 0; i < numCurves; ++i) {
    const auto &polys = orientedJordanCurves[i].getPolys();
    const auto &knots = orientedJordanCurves[i].getKnots();
    somePoints[i] = polys[0]((knots[0] + knots[1]) / 2);
    boxes[i] = boundingBox(orientedJordanCurves[i]);
    boundedness[i] = ::isBounded(orientedJordanCurves[i], 0);
  }
  for (int i = 0; i < numCurves - 1; ++i) {
    for (int j = i + 1; j < numCurves; ++j) {
      if (boxes[i].contain(boxes[j], tol)) candidates[i].push_back(j);
      if (boxes[j].contain(boxes[i], tol)) candidates[j].push_back(i);
    }
  }
  //
  PointsLocater locater(tol);
  for (int i = 0; i < numCurves; ++i) {
    std::vector<Vec<Real, 2>> queries;
    for (int k : candidates[i]) queries.push_back(somePoints[k]);
    auto loc = locater.operator()<Order>({orientedJordanCurves[i]}, queries,
                                         boundedness[i]);
    for (std::size_t j = 0; j < candidates[i].size(); ++j) {
      assert(loc[j] != 0);
      if (boundedness[i] == (loc[j] == 1)) {
        mat(i, candidates[i][j]) = 1;
        mat(candidates[i][j], i) = -1;
      }
    }
  }
  // step 2 : obtain the Hasse diagram from the inclusion matrix
  diagram.resize(numCurves + 1);
  diagram[numCurves].depth = -2;
  diagram[numCurves].parent = -1;
  std::vector<int> numAnc(numCurves);  // number of ancestors
  std::vector<int> parent(numCurves);
  for (int i = 0; i < numCurves; ++i) {
    numAnc[i] = (int)(std::count(&mat(0, i), &mat(0, i + 1), 1));
    if (numAnc[i] == 0) {
      parent[i] = numCurves;
      diagram[numCurves].depth = (boundedness[i]) ? (-1) : (0);
    }
  }
  assert(diagram[numCurves].depth != -2);
  // topological sort
  int numRemain = numCurves;
  int numPositive = 0;
  while (numRemain--) {
    int j =
        (int)(std::find(numAnc.cbegin(), numAnc.cend(), 0) - numAnc.cbegin());
    numAnc[j] = -1;  // so that j will never be selected again
    diagram[j].parent = parent[j];
    diagram[j].depth = diagram[parent[j]].depth + 1;
    if (diagram[j].depth % 2 == 0) ++numPositive;
    diagram[parent[j]].children.push_back(j);
    // handle the children
    for (int k = 0; k < numCurves; ++k) {
      if (mat(j, k) == 1) {
        parent[k] = j;
        --numAnc[k];
      }
    }
  }
  // step 3 : calculate the Betti numbers.
  if (diagram[numCurves].depth == -1) {
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
YinSet<2, 2> intersect(const YinSet<2, 2> &lhs, const YinSet<2, 2> &rhs,
                       Real tol) {
  auto segmented = meet(lhs, rhs, tol);
  return YinSet<2, 2>(segmented, tol);
}

//============================================================

bool equal(const Curve<2, 2> &lhs, const Curve<2, 2> &rhs, Real tol) {
  std::size_t N = lhs.getPolys().size();
  if (N != rhs.getPolys().size()) return false;
  VecCompare<Real, 2> vCmp(tol);
  auto findTopmostIdx = [&vCmp](const Curve<2, 2> &G) {
    Vec<Real, 2> r(-std::numeric_limits<Real>::max());
    std::size_t a;
    const auto &polys = G.getPolys();
    for (std::size_t i = 0; i < polys.size(); ++i)
      if (vCmp.compare(polys[i][0], r) == -1) {
        r = polys[i][0];
        a = i;
      }
    return a;
  };
  std::size_t k1 = findTopmostIdx(lhs);
  std::size_t k2 = findTopmostIdx(rhs);
  for (std::size_t i = 0; i < N; ++i) {
    if (vCmp(lhs.getPolys()[k1][0], rhs.getPolys()[k2][0]) != 0) return false;
    k1 = (k1 + 1) % N;
    k2 = (k2 + 1) % N;
  }
  return true;
}

template <>
bool YinSet<2, 2>::equal(const YinSet<2, 2> &rhs, Real tol) const {
  VecCompare<Real, 2> vCmp(tol);
  using Item = std::pair<int, bool>;
  std::map<Vec<Real, 2>, std::vector<Item>, VecCompare<Real, 2>> lookup(vCmp);

  auto findTopmost = [&vCmp](const Curve<2, 2> &G) {
    Vec<Real, 2> r(-std::numeric_limits<Real>::max());
    const auto &polys = G.getPolys();
    for (std::size_t i = 0; i < polys.size(); ++i)
      if (vCmp.compare(polys[i][0], r) == -1) r = polys[i][0];
    return r;
  };
  for (std::size_t i = 0; i < orientedJordanCurves.size(); ++i) {
    Vec<Real, 2> topmost = findTopmost(orientedJordanCurves[i]);
    auto ret = lookup.insert(std::make_pair(topmost, std::vector<Item>()));
    ret.first->second.push_back(std::make_pair(i, false));
  }
  for (const auto &j : rhs.orientedJordanCurves) {
    Vec<Real, 2> topmost = findTopmost(j);
    auto iter = lookup.find(topmost);
    if (iter == lookup.end()) return false;
    std::vector<Item> &items = iter->second;
    auto k = items.begin();
    for (; k != items.end(); ++k) {
      if (!k->second && ::equal(j, orientedJordanCurves[k->first], tol)) {
        k->second = true;
        break;
      }
    }
    if (k == items.end()) return false;
  }
  return true;
}

//============================================================

template <int Order>
std::string YinSet<2, Order>::getHasseString() const {
  if (orientedJordanCurves.empty()) return "YinSet empty.";
  const int w = 8;
  std::ostringstream oss;
  oss << std::left << std::setw(w) << " ";
  oss << std::left << std::setw(w) << "Parent";
  oss << std::left << std::setw(w) << "Depth";
  oss << std::left << std::setw(w) << "Orient";
  oss << std::left << "Children"
      << "\n";

  for (std::size_t i = 0; i < diagram.size(); ++i) {
    oss << std::left << std::setw(w) << i;
    oss << std::left << std::setw(w) << diagram[i].parent;
    oss << std::left << std::setw(w) << diagram[i].depth;
    if (i < orientedJordanCurves.size()) {
      oss << std::left << std::setw(w) << getOrientation(i);
    } else {
      oss << std::left << std::setw(w) << " ";
    }
    oss << "{";
    const std::vector<int> &children = diagram[i].children;
    for (int j : children) oss << j << ",";
    oss << "}"
        << "\n";
  }
  return oss.str();
}

template <int Order>
void YinSet<2, Order>::dump(std::ostream &os) const {
  const int N = orientedJordanCurves.size();
  os.write((char *)&N, sizeof(int));
  // save the boundaries according to the order of BFS
  const auto &rootOfForest = diagram.back();
  std::queue<int> Q;
  for (int i : rootOfForest.children) Q.push(i);
  while (!Q.empty()) {
    int k = Q.front();
    Q.pop();
    for (int i : diagram[k].children) Q.push(i);
    orientedJordanCurves[k].dump(os);
  }
}

template <int Order>
void YinSet<2, Order>::resetAllKinks(std::vector<Vertex> vertices) {
  sort(vertices.begin(), vertices.end());
  kinks = SimplicialComplex<Vertex>();
  auto numCurves = orientedJordanCurves.size();
  for (size_t i = 0; i < numCurves; ++i) {
    auto start =
        std::lower_bound(vertices.begin(), vertices.end(), Vertex(i, 0));
    auto end =
        std::lower_bound(vertices.begin(), vertices.end(), Vertex(i + 1, 0));
    for (auto Iter = start; Iter != end; ++Iter) {
      kinks.insert(Simplex<Vertex>{std::initializer_list<Vertex>{*Iter}});
    }
    reFitCurve(i);
  }
}

template <int Order>
int YinSet<2, Order>::insertKink(const Vertex &index) {
  auto sim = Simplex<Vertex>{std::initializer_list<Vertex>{index}};
  if (kinks.contain(sim)) return 0;
  kinks.insert(sim);
  reFitCurve(index.first);
  return 1;
}

template <int Order>
int YinSet<2, Order>::eraseKink(const Vertex &index) {
  auto sim = Simplex<Vertex>{std::initializer_list<Vertex>{index}};
  if (!kinks.contain(sim)) return 0;
  kinks.erase(sim);
  reFitCurve(index.first);
  return 1;
}

template <int Order>
void YinSet<2, Order>::reFitCurve(size_t i) {
  // Fit Curve is done in Mars.
  // if (kinks.getSimplexes().empty()) return;
  // const auto &sims = kinks.getSimplexes()[0];
  // SimplicialComplex<Vertex> tmp;
  // vector<Vec<Real, 2>> points;
  // auto start = sims.lower_bound(
  //          Simplex<Vertex>{std::initializer_list<Vertex>{{i, 0}}});
  // auto end = sims.lower_bound(
  //          Simplex<Vertex>{std::initializer_list<Vertex>{{i + 1, 0}}});
  // while (start != end) {
  //   tmp.insert(Simplex<Vertex>{
  //       std::initializer_list<Vertex>{{0,
  //       start->vertices.begin()->second}}});
  //   start++;
  // }
  // for (auto t : orientedJordanCurves[i].getKnots())
  //   points.push_back(orientedJordanCurves[i](t));
  // orientedJordanCurves[i].define(points, tmp);
}

template <int Order>
vector<Curve<2, Order>> YinSet<2, Order>::getSmoothCurves(Real tol) const {
  vector<Curve<2, Order>> res;
  const auto &sims = kinks.getSimplexes()[0];
  auto numJordanCurve = orientedJordanCurves.size();
  for (auto i = 0UL; i < numJordanCurve; ++i) {
    vector<Real> brks;
    auto start = sims.lower_bound(
        Simplex<Vertex>{std::initializer_list<Vertex>{{i, 0}}});
    auto end = sims.lower_bound(
        Simplex<Vertex>{std::initializer_list<Vertex>{{i + 1, 0}}});
    while (start != end) {
      const auto &index = *start->vertices.begin();
      brks.push_back(
          orientedJordanCurves[index.first].getKnots()[index.second]);
      ++start;
    }
    orientedJordanCurves[i].split(brks, res, tol);
  }
  return res;
}

template <int Order>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
auto YinSet<2, Order>::cutCell(const Box<Dim> &box, const Interval<Dim> &range,
                               bool addInner) const
    -> std::tuple<vector<vector<YinSetPtr>>,
                  vector<vector<vector<Curve<2, Order>>>>, vector<vector<int>>> {
  const int N0 = box.hi()[0] - box.lo()[0] + 1;
  const int N1 = box.hi()[1] - box.lo()[1] + 1;
  rVec size(box.size());
  vector<vector<YinSetPtr>> ret(N0, vector<YinSetPtr>(N1, nullptr));
  // since in pastCells will attach curve boundary to grid line, and will
  // increase tol error.
  Real tol = distTol() / 2;

  auto lo = range.lo();
  auto hi = range.hi();
  auto h = (hi - lo) / size;

  // calculate the intersections' parameters
  auto intersections = CutCellHelper<Order>::intersectGridLine(
      lo, hi, h, orientedJordanCurves, newtonTol());

  vector<vector<vector<Curve<2, Order>>>> gridCurves(
      N0, vector<vector<Curve<2, Order>>>(N1));
  // for (auto v : intersections[0]) {
  //   auto p = orientedJordanCurves[0](v);
  //   auto r = Vec<Real, 2>((p - lo) / h);
  //   std::cout << "p = " << p << ", r = " << r << std::endl;
  // }
  // partition the curves to vector<vector<Curve<Dim, Order>>>
  CutCellHelper<Order>::splitCurves(lo, h, intersections, orientedJordanCurves,
                                    gridCurves, newtonTol() *  100);

  // past in every cells.
  CutCellHelper<Order>::pastCells(lo, h, box, gridCurves, ret, tol);

  // fill inner cell rectangles.
  auto tags = CutCellHelper<Order>::fillInner(lo, h, box, *this, gridCurves,
                                              ret, tol, addInner);

  return {std::move(ret), std::move(gridCurves), std::move(tags)};
}

//============================================================
template class YinSet<2, 2>;
template class YinSet<2, 4>;
