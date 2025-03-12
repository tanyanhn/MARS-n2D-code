#include "Marsn2D/InterfaceGraph.h"

#include <deque>

#include "Core/VecCompare.h"

namespace MARSn2D {

InterfaceGraph::InterfaceGraph(
    vector<EdgeMark>&& edges,
    const vector<SmoothnessIndicator>& smoothConditions)
    : edges_(edges) {
  decomposeEdgeSet(smoothConditions, edges_, trials_, circuits_);
}
void InterfaceGraph::decomposeEdgeSet(
    const vector<SmoothnessIndicator>& smoothConditions,
    const vector<EdgeMark>& edges, vector<vector<EdgeIndex>>& trials,
    vector<vector<EdgeIndex>>& circuits) {
  VecCompare<Real, DIM> vCmp(distTol());
  unordered_map<int, int> prevSmoothPairs;
  unordered_map<int, int> postSmoothPairs;
  for (const auto& [pre, pos] : smoothConditions) {
    if (vCmp(edges[pre - 1].back(), edges[pos - 1].front()) != 0) {
      throw std::runtime_error("SmoothCondition is not continue.");
    }
    prevSmoothPairs[pre - 1] = pos - 1;
    postSmoothPairs[pos - 1] = pre - 1;
  }

  vector<int> edgeMarks(edges.size(), 0);
  for (int i = 0; i < edges.size(); ++i) {
    if (edgeMarks[i] == 0) {
      std::deque<int> line;
      line.push_back(i);
      edgeMarks[i] = 1;
      auto iter = prevSmoothPairs.find(line.back());
      while (iter != prevSmoothPairs.end() && edgeMarks[iter->second] == 0) {
        line.push_back(iter->second);
        edgeMarks[iter->second] = 1;
        iter = prevSmoothPairs.find(line.back());
      }
      iter = postSmoothPairs.find(line.front());
      while (iter != postSmoothPairs.end() && edgeMarks[iter->second] == 0) {
        line.push_front(iter->second);
        edgeMarks[iter->second] = 1;
        iter = postSmoothPairs.find(line.front());
      }

      if (iter != postSmoothPairs.end()) {
        if (iter->first == line.front() && iter->second == line.back())
          circuits.emplace_back(line.begin(), line.end());
        else
          throw std::runtime_error("Exist edge used multiple times.");
      } else {
        trials.emplace_back(line.begin(), line.end());
      }
    }
  }
}

template <int Order>
approxInterfaceGraph<Order>::approxInterfaceGraph(
    vector<EdgeMark>&& edgeMarks,
    const vector<SmoothnessIndicator>& smoothConditions,
    vector<vector<EdgeIndex>>&& cyclesEdgesId,
    vector<vector<size_t>>&& YinSetId, Real tol)
    : undirectGraph(std::move(edgeMarks), smoothConditions),
      cyclesEdgesId_(std::move(cyclesEdgesId)),
      yinSetId_(std::move(YinSetId)),
      tol_(tol) {
  updateCurve();
}

template <int Order>
void approxInterfaceGraph<Order>::updateCurve() {
  auto& marks_ = undirectGraph.edges_;
  auto& trials_ = undirectGraph.trials_;
  auto& circuits_ = undirectGraph.circuits_;
  edges_.resize(marks_.size());
  Real tol = this->tol_;

  auto fitFunction = [&marks_, tol](const vector<EdgeIndex>& spline,
                                    Curve<2, Order>::BCType type,
                                    vector<Edge>& output) {
    vector<Vertex> knots;
    vector<size_t> brksId;
    vector<Real> brks;
    VecCompare<Real, DIM> vCmp(tol);
    brksId.push_back(0);
    for (auto crvId : spline) {
      if (!knots.empty()) {
        auto p0 = knots.back();
        auto p1 = marks_[crvId].front();
        if (vCmp.compare(p0, p1) != 0)
          throw std::runtime_error("Invalid edge id");
        knots.pop_back();
        p0 = knots.back();
        if (vCmp.compare(p0, p1) == 0)
          throw std::runtime_error("Invalid edge id");
      }
      knots.insert(knots.end(), marks_[crvId].begin(), marks_[crvId].end());
      brksId.push_back(knots.size() - 1);
    }

    auto crv = fitCurve<Order>(knots, type);
    for (auto brkId : brksId) brks.push_back(crv.getKnots()[brkId]);

    vector<Edge> pieces;
    crv.split(brks, pieces, tol);

    for (size_t i = 0; i < pieces.size(); ++i)
      output[spline[i]] = std::move(pieces[i]);
  };
  for (const auto& trial : trials_)
    fitFunction(trial, Curve<DIM, Order>::notAKnot, edges_);

  for (const auto& circuit : circuits_)
    fitFunction(circuit, Curve<DIM, Order>::periodic, edges_);
}

template <int Order>
auto approxInterfaceGraph<Order>::approxJordanCurves() const
    -> vector<OrientedJordanCurve> {
  vector<OrientedJordanCurve> ret;

  for (const auto& cycle : cyclesEdgesId_) {
    Edge edge;
    for (const auto& edgeId : cycle) {
      if (edgeId > 0) {
        edge.concat(edges_[edgeId - 1]);
      } else if (edgeId < 0) {
        edge.concat(edges_[-edgeId - 1].reverse());
      } else {
        throw std::runtime_error("invalid edge id");
      }
    }
    ret.emplace_back(std::move(edge));
  }

  return ret;
}

template <int Order>
auto approxInterfaceGraph<Order>::approxYinSet() const
    -> vector<YinSet<DIM, Order>> {
  vector<YinSet<DIM, Order>> ret;

  auto jordanCurves = approxJordanCurves();

  for (const auto& vecId : yinSetId_) {
    std::vector<OrientedJordanCurve> boundary;
    boundary.reserve(vecId.size());
    for (auto id : vecId) {
      boundary.emplace_back(std::move(jordanCurves[id]));
    }
    ret.emplace_back(std::move(boundary), tol_);
  }

  return ret;
}

template <int Order>
auto approxInterfaceGraph<Order>::accessEdges()
    -> vector<std::pair<typename vector<Edge>::iterator,
                        typename vector<EdgeMark>::iterator>> {
  using std::vector;
  vector<std::pair<typename vector<Edge>::iterator,
                   typename vector<EdgeMark>::iterator>>
      ret;

  auto& marks = undirectGraph.edges_;
  auto iterMarks = marks.begin();
  auto iterEdges = edges_.begin();
  for (; iterMarks != marks.cend(); ++iterMarks, ++iterEdges) {
    ret.emplace_back(iterEdges, iterMarks);
  }

  return ret;
}

//===============================================

template class approxInterfaceGraph<2>;
template class approxInterfaceGraph<4>;

}  // namespace MARSn2D