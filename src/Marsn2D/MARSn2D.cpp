#include "Marsn2D/MARSn2D.h"

#include "Recorder/ProgressBar.h"
#include "Recorder/Timer.h"

namespace Marsn2D {

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::discreteFlowMap(const VectorFunction<2> &v,
                                                    EdgeMark &marks, Real tn,
                                                    Real dt) const {
  Timer t("discreteFlowMap");
  // using Eigen.
  if (dt > 0)
    TI_->timeStep(v, marks, tn, dt);
  // for (auto & pt : marks) {
  //   pt = TI_->timeStep(v, pt, tn, dt);
  // }
}

template <int Order, template <int> class VelocityField>
OPTNONE_FUNC void MARSn2D<Order, VelocityField>::plot(
    const IG &ig, int step, const PlotConfig &plotConfig) const {
  auto yinSets = ig.approxYinSet();
  std::string fileName = plotConfig.fName + "_Step" + std::to_string(step);
  std::string tail = (plotConfig.output == NORMAL) ? "_f.dat" : "_c.dat";
  std::ofstream of(fileName + tail, std::ios_base::binary);
  int num = yinSets.size();
  if (plotConfig.output == NORMAL) {
    of.write((char *)&num, sizeof(num));
  }
  for (int i = 0; i <= num - 1; i++) {
    auto &yinSet = yinSets[i];
    if (plotConfig.output == NORMAL) {
      yinSet.dump(of);
    } else if (plotConfig.output == CUTCELL) {
      auto [res, boundary, tags] = yinSet.cutCell(
          plotConfig.box, plotConfig.range, plotConfig.plotInner);
      dumpVecYinSet<Order>(res, of);
    }
  }
}

template <int Order, template <int> class VelocityField>
OPTNONE_FUNC void MARSn2D<Order, VelocityField>::trackInterface(
    const VelocityField<DIM> &v, IG &ig, Real StartTime, Real dt, Real EndTime,
    const PlotConfig &plotConfig) const {
  Real tn = StartTime;
  int stages = 0;
  int polyStep = 0;
  Real k = 0;
  if (StartTime > EndTime) {
    stages = ceil(abs(EndTime - StartTime) / abs(dt));
    k = (EndTime - StartTime) / stages;
    polyStep = stages / plotConfig.opStride;
  }
  int step = 0;
  vector<vector<int>> marksHistory;
  vector<vector<Real>> lengthHistory;
  // marksHistory.push_back(ig.countMarks());
  // lengthHistory.push_back(ig.countLengths());
  ProgressBar bar(stages, "Tracking...");
  timeStep(v, ig, tn, 0);
  while (step < stages) {
    // if (plotConfig.output != NONE && (step + 40) % polyStep == 0) {
    //   plot(ig, step, plotConfig);
    // }
    if (plotConfig.output != NONE && step % polyStep == 0) {
      plot(ig, step, plotConfig);
    }
    if (printDetail)
      std::cout << "Step: " << step << "     time now: " << tn << '\n';
    else
      bar.update();
    if (step == stages - 1) {
      timeStep(v, ig, tn, EndTime - tn);
      tn = EndTime;
    } else {
      timeStep(v, ig, tn, k);
      tn += k;
    }
    std::cout.flush();
    marksHistory.push_back(ig.countMarks());
    lengthHistory.push_back(ig.countLengths());
    step++;
  }
  if (plotConfig.output != NONE) {
    plot(ig, step, plotConfig);
    auto dir = getExportDir();
    std::ofstream of1(plotConfig.fName + "_00marksHistory" + ".dat",
                      std::ios_base::binary);
    std::ofstream of2(plotConfig.fName + "_00LengthHistory" + ".dat",
                      std::ios_base::binary);
    int rows = (int)marksHistory[0].size();
    int cols = (int)marksHistory.size();
    of1.write((char *)&rows, sizeof(int));
    of1.write((char *)&cols, sizeof(int));
    of2.write((char *)&rows, sizeof(int));
    of2.write((char *)&cols, sizeof(int));
    for (size_t i = 0; i < marksHistory.size(); i++) {
      of1.write((char *)marksHistory[i].data(), rows * sizeof(int));
      of2.write((char *)lengthHistory[i].data(), rows * sizeof(Real));
    }
  }
}

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::locateLongEdges(
    EdgeMark &marks, const vector<Real> &hL,
    vector<std::pair<unsigned int, unsigned int>> &indices2Num,
    Real efficientOfHL) const {
  size_t numEdge = marks.size() - 1;
  for (size_t i = 0; i < numEdge; i++) {
    Real length = norm(marks[i + 1] - marks[i]);
    // insert when both side curv not satisfy condition.
    Real constrain = efficientOfHL * std::max(hL[i], hL[i + 1]);
    if (length > constrain) {
      unsigned int num = (int)std::ceil(length / constrain) * 2;
      indices2Num.emplace_back(i, num);
    }
  }
}

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::locateTinyEdges(
    EdgeMark &marks, const vector<Real> &hL, vector<unsigned int> &indices,
    Real efficientOfHL) const {
  size_t numEdge = marks.size() - 1;
  vector<Real> lengths;
  for (size_t i = 0; i < numEdge; i++) {
    Real length = norm(marks[i + 1] - marks[i]);
    lengths.push_back(length);
  }
  for (size_t i = 0; i < numEdge; i++) {
    // remove when both side curv not satisfy condition.
    Real constrain = efficientOfHL * std::min(hL[i], hL[i + 1]);
    if (lengths[i] < constrain) {
      if (i == 0) {
        indices.emplace_back(++i);
      } else if (i == numEdge - 1) {
        indices.emplace_back(i);
      } else {
        if (lengths[i - 1] / std::min(hL[i], hL[i - 1]) <
            lengths[i + 1] / std::min(hL[i + 1], hL[i + 2])) {
          indices.emplace_back(i);
        } else {
          indices.emplace_back(++i);
        }
      }
    }
  }
}

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::insertMarks(
    const VelocityField<2> &v, Real preT, Real dt,
    const vector<std::pair<unsigned int, unsigned int>> &indices2Num,
    EdgeMark &marks, Edge &preEdge, vector<Real> &curv,
    vector<Real> &para) const {
  const auto &polys = preEdge.getPolys();
  const auto &knots = preEdge.getKnots();
  EdgeMark newMarks;
  vector<Polynomial<Order, Point>> newPolys;
  vector<Real> newKnots;
  vector<Real> newCurv;
  newMarks.reserve(marks.size() + indices2Num.size());
  newPolys.reserve(polys.size() + indices2Num.size());
  newKnots.reserve(knots.size() + indices2Num.size());
  if (curvConfig_.used) newCurv.reserve(curv.size() + indices2Num.size());
  size_t numEdge = marks.size() - 1;
  size_t index = 0;
  for (size_t i = 0; i < numEdge; i++) {
    newMarks.push_back(marks[i]);
    newPolys.push_back(polys[i]);
    newKnots.push_back(knots[i]);
    if (curvConfig_.used) newCurv.push_back(curv[i]);
    if (index < indices2Num.size()) {
      auto [j, num] = indices2Num[index];
      if (i == j) {
        index++;
        EdgeMark subMarks;
        const auto &poly = polys[i];
        const Real s = knots[i + 1] - knots[i];
        const Real ds = s / num;
        for (size_t k = 1; k < num; k++) {
          Real t = (Real)k * ds;
          Point pt = poly(t);
          subMarks.push_back(pt);
          newPolys.push_back(poly.translate(t));
          newKnots.push_back(knots[i] + t);
          if (curvConfig_.used)
            newCurv.push_back(Edge::curvature(getComp(newPolys.back(), 0),
                                              getComp(newPolys.back(), 1), 0));
        }
        discreteFlowMap(v, subMarks, preT, dt);
        newMarks.insert(newMarks.end(), subMarks.begin(), subMarks.end());
      }
    }
  }
  newMarks.push_back(marks.back());
  Real t = knots.back() - newKnots.back();
  newKnots.push_back(knots.back());
  if (curvConfig_.used)
    newCurv.push_back(Edge::curvature(getComp(newPolys.back(), 0),
                                      getComp(newPolys.back(), 1), t));
  marks = std::move(newMarks);
  para = newKnots;
  preEdge = std::move(Edge(std::move(newKnots), std::move(newPolys)));
  if (curvConfig_.used) curv = std::move(newCurv);
}

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::removeMarks(
    const vector<unsigned int> &indices, EdgeMark &marks, vector<Real> &curv,
    vector<Real> &para) const {
  EdgeMark newMarks;
  vector<Real> newCurv;
  vector<Real> newPara;
  newMarks.reserve(marks.size() - indices.size());
  if (curvConfig_.used) newCurv.reserve(curv.size() - indices.size());
  size_t numMarks = marks.size();
  size_t index = 0;
  bool continueFlag = false;
  for (size_t i = 0; i < numMarks; i++) {
    if (index == indices.size() || i != indices[index]) {
      newMarks.push_back(marks[i]);
      newPara.push_back(para[i]);
      if (curvConfig_.used) newCurv.push_back(curv[i]);
      continueFlag = false;
    } else {
      if (continueFlag) {
        newMarks.push_back(marks[i]);
        newPara.push_back(para[i]);
        if (curvConfig_.used) newCurv.push_back(curv[i]);
        continueFlag = false;
      } else {
        continueFlag = true;
      }
      ++index;
    }
  }
  marks = std::move(newMarks);
  para = std::move(newPara);
  if (curvConfig_.used) curv = std::move(newCurv);
}

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::updateHL(const vector<Real> &curv,
                                             const EdgeMark &marks,
                                             vector<Real> &hL) const {
  if (!curvConfig_.used) {
    hL.resize(marks.size(), hL_);
    return;
  }
  hL.resize(curv.size());
  Real curvMax = curv.front();
  Real curvMin = curv.front();
  for (const auto &c : curv) {
    curvMax = std::max(curvMax, c);
    curvMin = std::min(curvMin, c);
  }
  curvMax = std::min(curvMax, curvConfig_.rhoCRange.hi()[0]);
  curvMin = std::max(curvMin, curvConfig_.rhoCRange.lo()[0]);
  Real rMin = std::max(curvConfig_.rCMin, curvMin / curvMax);

  for (size_t i = 0; i < curv.size(); i++) {
    if (curv[i] <= curvMin)
      hL[i] = hL_;
    else if (curv[i] >= curvMax)
      hL[i] = hL_ * rMin;
    else {
      hL[i] = rMin * hL_ + (1 - rMin) * hL_ *
                               curvConfig_.sigma((1 / curv[i] - 1 / curvMax) /
                                                 (1 / curvMin - 1 / curvMax));
    }
  }
}

template <int Order, template <int> class VelocityField>
OPTNONE_FUNC void MARSn2D<Order, VelocityField>::averageSplit(
    const VelocityField<2> &v, Real preT, Real dt, EdgeMark &marks,
    const Edge &preEdge, vector<Real> &hL, vector<Real> &para,
    Real efficientOfHL, int notaKnotLocation) const {
  Timer timer("averageSplit()");
  //TODO(ytan): need more split on pp_1 to prevent >hL edge.
  auto locateSplit = [&v, &preEdge, preT, dt, this](auto &p0, auto &p1, auto hL,
                                                    auto t0,
                                                    auto t1) OPTNONE_FUNC {
    Real localT0 = t0;
    Real localT1 = t1;
    auto midT = (localT0 + localT1) / 2;
    auto midP = EdgeMark{preEdge(midT)};
    discreteFlowMap(v, midP, preT, dt);
    int whileCount = 20;
    Real d0 = norm(midP[0] - p0);
    Real d1 = norm(midP[0] - p1);
    Real tmp;
    while (d0 > d1 || d0 < hL) {
      tmp = midT;
      if (d0 > d1) {
        midT = (localT0 + midT) / 2;
        localT1 = tmp;
      } else if (d0 < hL) {
        midT = (midT + localT1) / 2;
        localT0 = tmp;
      } else {
        throw std::runtime_error("averageSplit unable find mid point.");
      }
      midP[0] = preEdge(midT);
      discreteFlowMap(v, midP, preT, dt);
      d0 = norm(midP[0] - p0);
      d1 = norm(midP[0] - p1);
      if (--whileCount <= 0) {
        throw std::runtime_error("averageSplit overloop.");
      }
    }
    return midP[0];
  };
  std::unordered_map<int, Vertex> splitVertexes;
  if (notaKnotLocation &
      approxInterfaceGraph<Order>::notAKnotBoundaryLocation::Left) {
    auto p0 = marks[0];
    auto p1 = marks[1];
    auto p2 = marks[2];
    if (norm(p1 - p0) > norm(p2 - p1)) {
      Real localHL = efficientOfHL * std::min(hL[0], hL[1]);
      // try {
        auto splitVertex = locateSplit(p0, p1, localHL, para[0], para[1]);
        splitVertexes.insert({-1, Vertex(splitVertex)});
      // } catch (std::runtime_error &e) {
      //   auto splitVertex = locateSplit(p0, p1, localHL, para[0], para[1]);
      // }
    }
  }
  if (notaKnotLocation &
      approxInterfaceGraph<Order>::notAKnotBoundaryLocation::Right) {
    auto p0 = marks[marks.size() - 1];
    auto p1 = marks[marks.size() - 2];
    auto p2 = marks[marks.size() - 3];
    if (norm(p1 - p0) > norm(p2 - p1)) {
      Real localHL =
          efficientOfHL * std::min(hL[marks.size() - 1], hL[marks.size() - 2]);
      // try {
        auto splitVertex = locateSplit(p0, p1, localHL, para[marks.size() - 1],
                                       para[marks.size() - 2]);
        splitVertexes.insert({1, Vertex(splitVertex)});
      // } catch (std::runtime_error &e) {
      //   auto splitVertex = locateSplit(p0, p1, localHL, para[marks.size() - 1],
      //                                  para[marks.size() - 2]);
      // }
    }
  }
  if (splitVertexes.empty()) return;
  vector<Vertex> newMarks;
  vector<Real> newHL;
  newMarks.push_back(marks[0]);
  newHL.push_back(hL[0]);
  if (splitVertexes.contains(-1)) {
    newMarks.push_back(splitVertexes.at(-1));
    newHL.push_back(std::min(hL[0], hL[1]));
  }
  newMarks.insert(newMarks.end(), marks.begin() + 1, marks.end() - 1);
  newHL.insert(newHL.end(), hL.begin() + 1, hL.end() - 1);
  if (splitVertexes.contains(1)) {
    newMarks.push_back(splitVertexes.at(1));
    newHL.push_back(std::min(hL[hL.size() - 2], hL.back()));
  }
  newMarks.push_back(marks.back());
  newHL.push_back(hL.back());
  marks = std::move(newMarks);
  hL = std::move(newHL);
}

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::timeStep(const VelocityField<DIM> &v,
                                             IG &ig, Real tn, Real dt) const {
  Timer t("timeStep");
  int insertCount = 0;
  int removeCount = 0;
  auto stepCrv = [this, &v, tn, dt, &insertCount, &removeCount](
                     Edge &preEdge, EdgeMark &marks,
                     int notaKnotLocation) OPTNONE_FUNC {
    vector<Real> curv;
    vector<Real> para = preEdge.getKnots();
    vector<Real> hL;
    Real ARMSDist = 1 * rTiny_;
    Real splitDist = ARMSDist * 1/3.0;
    discreteFlowMap(v, marks, tn, dt);

    if (curvConfig_.used) {
      curv.resize(marks.size());
      const auto &polys = preEdge.getPolys();
      const auto &knots = preEdge.getKnots();
      for (size_t i = 0; i < marks.size() - 1; i++) {
        curv[i] =
            Edge::curvature(getComp(polys[i], 0), getComp(polys[i], 1), 0);
      }
      curv.back() =
          Edge::curvature(getComp(polys.back(), 0), getComp(polys.back(), 1),
                          knots[knots.size() - 1] - knots[knots.size() - 2]);
    }

    bool inserted = true;
    while (inserted && insertCount < insertTimesMax) {
      updateHL(curv, marks, hL);
      vector<std::pair<unsigned int, unsigned int>> indices2Num;
      /** TODO(ytan) need efficientOfHL = 1 - 2 * rTiny_ fill hL upper bound,
       *  since one can delete | rTiny | 1 | rTiny|rightend make 1 + 2 *
       *rTiny_ edge.
       **/
      // locateLongEdges(marks, hL, indices2Num, 1 - 2 * ARMSDist);
      locateLongEdges(marks, hL, indices2Num, 1 - 2 * ARMSDist);
      inserted = !indices2Num.empty();
      if (inserted) {
        insertMarks(v, tn, dt, indices2Num, marks, preEdge, curv, para);
        insertCount++;
        if (insertCount >= 10) 
          std::cout << insertCount << ' ' << marks.size() << std::endl;
        if (insertCount >= insertTimesMax)
          throw std::runtime_error("insertCount >= insertTimesMax");
      }
    }

    bool removed = true;
    while (removed && removeCount < removeTimesMax && marks.size() > 2) {
      updateHL(curv, marks, hL);
      vector<unsigned int> indices;
      locateTinyEdges(marks, hL, indices, ARMSDist);
      removed = !indices.empty();
      if (removed) {
        removeMarks(indices, marks, curv, para);
        removeCount++;
        if (removeCount >= removeTimesMax)
          throw std::runtime_error("removeCount >= removeTimesMax");
      }
    }

    // split boundary edge for not-a-knot splines.
    updateHL(curv, marks, hL);
    averageSplit(v, tn, dt, marks, preEdge, hL, para, splitDist,
                 notaKnotLocation);

#ifndef NDEBUG
    if (!checkMarks(marks, hL)) {
      std::cout << "Marks are not valid.\n";
      exit(1);
    }
#endif  // !NDEBUG
  };

  auto mark_edge = ig.accessEdges();
  // #pragma omp parallel for default(shared) schedule(static)
  for (auto [edgeIter, markIter, notaKnotLocation] : mark_edge) {
    insertCount = 0;
    removeCount = 0;
    // discreteFlowMap(v, *markIter, tn, dt);
    stepCrv(*edgeIter, *markIter, notaKnotLocation);
    if (printDetail) {
      std::cout << fmt::v11::format("Insert: {}, Remove: {}. \n", insertCount,
                                    removeCount);
    }
  }
  // #pragma omp critical
  ig.updateCurve();
}

template <int Order, template <int> class VelocityField>
bool MARSn2D<Order, VelocityField>::checkMarks(const EdgeMark &marks,
                                               const vector<Real> &hL) const {
  size_t numEdge = marks.size() - 1;
  for (size_t i = 0; i < numEdge; i++) {
    Real length = norm(marks[i + 1] - marks[i]);
    Real upBound = std::max(hL[i], hL[i + 1]);
    Real loBound = std::min(hL[i], hL[i + 1]);
    if (length < rTiny_ * loBound) {
      std::cout << "Edge " << i << " is too short.\n";
      return false;
    }
    if (length > upBound) {
      std::cout << "Edge " << i << " is too long.\n";
      return false;
    }
  }
  return true;
}

//===========================================

template class MARSn2D<2, VectorFunction>;
template class MARSn2D<4, VectorFunction>;

}  // namespace Marsn2D