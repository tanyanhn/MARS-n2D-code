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
  // TI_->timeStep(v, marks, tn, dt);
  for (auto & pt : marks) {
    pt = TI_->timeStep(v, pt, tn, dt);
  }
}

template <int Order, template <int> class VelocityField>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
void MARSn2D<Order, VelocityField>::plot(const IG &ig, int step,
                                         const PlotConfig &plotConfig) const {
  auto yinSets = ig.approxYinSet();
  std::string fileName = plotConfig.fName + "_Step" + std::to_string(step);
  std::ofstream of(fileName + "_c.dat", std::ios_base::binary);
  int num = yinSets.size();
  for (int i = num - 1; i >= 0; i--) {
    auto& yinSet = yinSets[i];
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
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
void MARSn2D<Order, VelocityField>::trackInterface(
    const VelocityField<DIM> &v, IG &ig, Real StartTime, Real dt, Real EndTime,
    const PlotConfig& plotConfig) const {
  Real tn = StartTime;
  int stages = ceil(abs(EndTime - StartTime) / abs(dt));
  Real k = (EndTime - StartTime) / stages;
  int step = 0;
  int polyStep = stages / plotConfig.opStride;
  vector<vector<int>> marksHistory;
  vector<vector<Real>> lengthHistory;
  marksHistory.push_back(ig.countMarks());
  lengthHistory.push_back(ig.countLengths());
  ProgressBar bar(stages, "Tracking...");
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
    std::ofstream of1(dir + "00marksHistory.dat", std::ios_base::binary);
    std::ofstream of2(dir + "00LengthHistory.dat", std::ios_base::binary);
    int rows = marksHistory[0].size();
    int cols = marksHistory.size();
    of1.write((char*)&rows, sizeof(int));
    of1.write((char*)&cols, sizeof(int));
    of2.write((char*)&rows, sizeof(int));
    of2.write((char*)&cols, sizeof(int));
    for (size_t i = 0; i < marksHistory.size(); i++) {
      of1.write((char*)marksHistory[i].data(), rows * sizeof(int));
      of2.write((char*)lengthHistory[i].data(), rows * sizeof(Real));
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
      unsigned int num = ceil(length / constrain) * 2;
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
    EdgeMark &marks, Edge &preEdge, vector<Real> &curv) const {
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
  preEdge = std::move(Edge(std::move(newKnots), std::move(newPolys)));
  if (curvConfig_.used) curv = std::move(newCurv);
}

template <int Order, template <int> class VelocityField>
void MARSn2D<Order, VelocityField>::removeMarks(
    const vector<unsigned int> &indices, EdgeMark &marks,
    vector<Real> &curv) const {
  EdgeMark newMarks;
  vector<Real> newCurv;
  newMarks.reserve(marks.size() - indices.size());
  if (curvConfig_.used) newCurv.reserve(curv.size() - indices.size());
  size_t numMarks = marks.size();
  size_t index = 0;
  bool continueFlag = false;
  for (size_t i = 0; i < numMarks; i++) {
    if (index == indices.size() || i != indices[index]) {
      newMarks.push_back(marks[i]);
      if (curvConfig_.used) newCurv.push_back(curv[i]);
      continueFlag = false;
    } else {
      if (continueFlag) {
        newMarks.push_back(marks[i]);
        if (curvConfig_.used) newCurv.push_back(curv[i]);
        continueFlag = false;
      } else {
        continueFlag = true;
      }
      ++index;
    }
  }
  marks = std::move(newMarks);
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
void MARSn2D<Order, VelocityField>::timeStep(const VelocityField<DIM> &v,
                                             IG &ig, Real tn, Real dt) const {
  Timer t("timeStep");
  int insertCount = 0;
  int removeCount = 0;
  auto stepCrv = [ this, &v, tn, dt, &insertCount, &
                   removeCount ](Edge & preEdge, EdgeMark & marks)
#ifdef OPTNONE
      __attribute__((optnone))
#endif  // OPTNONE
  {
    vector<Real> curv;
    vector<Real> hL;
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
       *  since one can delete | rTiny | 1 | rTiny|rightend make 1 + 2 * rTiny_
       *  edge.
       **/
      locateLongEdges(marks, hL, indices2Num, 1 - rTiny_);
      inserted = !indices2Num.empty();
      if (inserted) {
        insertMarks(v, tn, dt, indices2Num, marks, preEdge, curv);
        insertCount++;
        if (insertCount >= insertTimesMax)
          throw std::runtime_error("insertCount >= insertTimesMax");
      }
    }

    bool removed = true;
    while (removed && removeCount < removeTimesMax && marks.size() > 2) {
      updateHL(curv, marks, hL);
      vector<unsigned int> indices;
      locateTinyEdges(marks, hL, indices, rTiny_);
      removed = !indices.empty();
      if (removed) {
        // if (indices.back() == marks.size() - 2) {
        //   indices.pop_back();
        //   if (indices.empty() || indices.back() != marks.size() - 3)
        //     indices.push_back(marks.size() - 3);
        // }
        removeMarks(indices, marks, curv);
        removeCount++;
        if (removeCount >= removeTimesMax)
          throw std::runtime_error("removeCount >= removeTimesMax");
      }
    }

#ifndef NDEBUG
    if (!checkMarks(marks, hL)) {
      std::cout << "Marks are not valid.\n";
      exit(1);
    }
#endif  // !NDEBUG
  };

  auto mark_edge = ig.accessEdges();
  // #pragma omp parallel for default(shared) schedule(static)
  for (auto [edgeIter, markIter] : mark_edge) {
    insertCount = 0;
    removeCount = 0;
    stepCrv(*edgeIter, *markIter);
    if (printDetail) {
      std::cout << std::format("Insert: {}, Remove: {}. \n", insertCount,
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