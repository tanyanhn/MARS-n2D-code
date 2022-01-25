#include "OrientedJordanCurve.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <istream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include "Core/Curve.h"
#include "SimplicialComplex.h"
#include "YinSet.h"
#include "YinSet/SegmentedRealizableSpadjor.h"

using std::string;

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(const string& parameters) {
  SimplicialComplex kinks;
  define(parameters, kinks);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(const string& parameters,
                                             SimplicialComplex& kinks) {
  assert(kinks.getNSim() == -1);
  std::stringstream iss(parameters);
  define(iss, kinks);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(std::istream& iss,
                                             SimplicialComplex& kinks) {
  unsigned int nPolys;
  // iss.read((char*)&nPolys, sizeof(unsigned int));
  iss >> nPolys;
  std::vector<Vec<Real, Dim>> points(nPolys + 1);
  for (size_t i = 0; i <= nPolys; i++) {
    // Real buf[Dim];
    // Real* pbuf = buf;
    // iss.read((char*)buf, sizeof(Real) * Dim);
    for (int j = 0; j < Dim; j++)
      // points[i][j] = *pbuf++;
      iss >> points[i][j];
  }
  unsigned int id;
  while (iss >> id)
    kinks.insert(Simplex{std::initializer_list<unsigned int>{id}});

  define(points, kinks);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(
    const std::vector<Vec<Real, Dim>>& points,
    const SimplicialComplex& kinks) {
  if (kinks.getSimplexes().size() == 0) {
    Curve<Dim, Order> res = fitCurve<Order>(points, periodic);
    this->knots = res.getKnots();
    this->polys = res.getPolys();
  } else {
    int pre = -1, pos = -1, last = -1;
    for (auto& simplex : kinks.getSimplexes()[0]) {
      if (pre == -1) {
        pre = *simplex.vertices.begin();
        last = pre;
      } else {
        pos = *simplex.vertices.begin();
        vector<Vec<Real, Dim>> subPoints(std::next(points.begin(), pre),
                                         std::next(points.begin(), pos + 1));
        Curve<Dim, Order> res = fitCurve<Order>(subPoints, nature);
        Real startT;
        if (this->knots.empty()) {
          startT = 0;
        } else {
          startT = this->knots.back();
          this->knots.pop_back();
        }
        for (auto t : res.getKnots()) {
          this->knots.push_back(startT + t);
        }
        this->polys.insert(this->polys.end(), res.getPolys().begin(),
                           res.getPolys().end());
        pre = pos;
      }
    }
    vector<Vec<Real, Dim>> subPoints(std::next(points.begin(), pre),
                                     points.end());
    if (last != 0)
      subPoints.insert(subPoints.end(), std::next(points.begin(), 1),
                       std::next(points.begin(), last + 1));
    Curve<Dim, Order> res = fitCurve<Order>(subPoints, nature);
    Real startT;
    if (this->knots.empty()) {
      startT = 0;
    } else {
      startT = this->knots.back();
      this->knots.pop_back();
    }
    for (auto t : res.getKnots()) {
      this->knots.push_back(startT + t);
    }
    this->polys.insert(this->polys.end(), res.getPolys().begin(),
                       res.getPolys().end());
  }
  return;
}

template <int Order>
void Circle<Order>::define(const std::string& parameters) {
  SimplicialComplex kinks;
  define(parameters, kinks);
}

template <int Order>
void Circle<Order>::define(const std::string& parameters,
                           SimplicialComplex& kinks) {
  assert(kinks.getNSim() == -1);
  std::stringstream iss(parameters);
  define(iss, kinks);
}

template <int Order>
void Circle<Order>::define(std::istream& iss, SimplicialComplex& kinks) {
  Vec<Real, 2> center;
  Real radius;
  bool orientation;
  Real hL;
  // iss.read((char*)&(center[0]), sizeof(Real));
  // iss.read((char*)&(center[1]), sizeof(Real));
  // iss.read((char*)&(radius), sizeof(Real));
  // iss.read((char*)&(orientation), sizeof(bool));
  // iss.read((char*)&(hL), sizeof(Real));
  iss >> center[0] >> center[1] >> radius >> orientation >> hL;
  int sign = orientation ? 1 : -1;
  size_t nPolys = radius * 2 * M_PI / hL + 1;
  Real dtheta = 2 * M_PI / nPolys;
  std::vector<Vec<Real, 2>> points(nPolys + 1);
  for (size_t i = 0; i <= nPolys; ++i) {
    points[i][0] = center[0] + radius * std::cos(sign * i * dtheta);
    points[i][1] = center[1] + radius * std::sin(sign * i * dtheta);
  }
  OrientedJordanCurve<2, Order>::define(points, kinks);
}

template <int Order>
void Rectangle<Order>::define(const std::string& parameters) {
  SimplicialComplex kinks;
  define(parameters, kinks);
}

template <int Order>
void Rectangle<Order>::define(const std::string& parameters,
                              SimplicialComplex& kinks) {
  assert(kinks.getNSim() == -1);
  std::stringstream iss(parameters);
  define(iss, kinks);
}

template <int Order>
void Rectangle<Order>::define(std::istream& iss, SimplicialComplex& kinks) {
  Vec<Real, 2> smallEnd, bigEnd;
  Real theta;
  bool orientation;
  Real hL;
  // iss.read((char*)&(smallEnd[0]), sizeof(Real));
  // iss.read((char*)&(smallEnd[1]), sizeof(Real));
  // iss.read((char*)&(bigEnd[0]), sizeof(Real));
  // iss.read((char*)&(bigEnd[1]), sizeof(Real));
  // iss.read((char*)&(theta), sizeof(Real));
  // iss.read((char*)&(orientation), sizeof(bool));
  // iss.read((char*)&(hL), sizeof(Real));
  iss >> smallEnd[0] >> smallEnd[1] >> bigEnd[0] >> bigEnd[1] >> theta >>
      orientation >> hL;
  size_t rowPolys = (bigEnd[0] - smallEnd[0]) / hL + 1,
         colPolys = (bigEnd[1] - smallEnd[1]) / hL + 1;
  Real dx = (bigEnd[0] - smallEnd[0]) / rowPolys,
       dy = (bigEnd[1] - smallEnd[1]) / colPolys;

  std::vector<Vec<Real, 2>> points;
  kinks.insert(Simplex{std::initializer_list<unsigned long>{points.size()}});
  for (size_t i = 0; i < rowPolys; ++i)
    points.push_back(Vec<Real, 2>{smallEnd[0] + i * dx, smallEnd[1]});
  kinks.insert(Simplex{std::initializer_list<unsigned long>{points.size()}});
  for (size_t i = 0; i < colPolys; ++i)
    points.push_back(Vec<Real, 2>{bigEnd[0], smallEnd[1] + i * dy});
  kinks.insert(Simplex{std::initializer_list<unsigned long>{points.size()}});
  for (size_t i = 0; i < rowPolys; ++i)
    points.push_back(Vec<Real, 2>{bigEnd[0] - i * dx, bigEnd[1]});
  kinks.insert(Simplex{std::initializer_list<unsigned long>{points.size()}});
  for (size_t i = 0; i < colPolys; ++i)
    points.push_back(Vec<Real, 2>{smallEnd[0], bigEnd[1] - i * dy});
  points.push_back(smallEnd);

  Real c = cos(theta), s = sin(theta), x, y;
  for (auto& p : points) {
    x = c * p[0] - s * p[1];
    y = s * p[0] + c * p[1];
    p[0] = x, p[1] = y;
  }

  OrientedJordanCurve<2, Order>::define(points, kinks);
}

template <int Order>
std::unique_ptr<OrientedJordanCurve<2, Order>>
CurveFactory<2, Order>::createCurve(const std::string& parameters) {
  SimplicialComplex kinks;
  return createCurve(parameters, kinks);
}

template <int Order>
std::unique_ptr<OrientedJordanCurve<2, Order>>
CurveFactory<2, Order>::createCurve(const std::string& parameters,
                                    SimplicialComplex& kinks) {
  assert(kinks.getNSim() == -1);
  std::stringstream iss(parameters);
  string type;
  iss >> type;
  std::unique_ptr<OrientedJordanCurve<2, Order>> ret(nullptr);
  if (type == "OrientedJordanCurve") {
    ret.reset(new OrientedJordanCurve<2, Order>());
    ret->define(iss, kinks);
  } else if (type == "Circle") {
    ret.reset(new Circle<Order>());
    ret->define(iss, kinks);
  } else if (type == "Rectangle") {
    ret.reset(new Rectangle<Order>());
    ret->define(iss, kinks);
  } else {
    assert(false && "undefined type.");
  }
  return ret;
}

template <int Order>
YinSet<2, Order> CurveFactory<2, Order>::createYinSet(
    const std::vector<std::string>& parameters) {
  size_t nCurves = parameters.size() - 1;
  std::vector<OrientedJordanCurve<2, Order>> curves(nCurves);
  Real tol = stod(parameters[0]);
  SimplicialComplex kinks;
  unordered_map<unsigned int, std::pair<unsigned int, unsigned int>>
      mVertex2Point;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> mPoint2Vertex;
  unsigned int id = 0;
  for (size_t i = 0; i < nCurves; ++i) {
    SimplicialComplex tmp;
    curves[i] = *createCurve(parameters[i + 1], tmp);
    if (tmp.getNSim() == 0) {
      for (auto& simplex : tmp.getSimplexes()[0]) {
        unsigned int j = *simplex.vertices.begin();
        auto index = std::make_pair(i, j);
        kinks.insert(Simplex{std::initializer_list<unsigned int>{id}});
        mVertex2Point[id] = index;
        mPoint2Vertex[index] = id;
        ++id;
      }
    }
  }
  SegmentedRealizableSpadjor<Order> segmentedSpadjor(curves, tol);

  YinSet<2, Order> ret(segmentedSpadjor, tol);
  ret.setSimplexes(kinks, mVertex2Point, mPoint2Vertex);

  return ret;
}

template class OrientedJordanCurve<2, 2>;
template class OrientedJordanCurve<2, 4>;
template struct Circle<2>;
template struct Circle<4>;
template struct Rectangle<2>;
template struct Rectangle<4>;
template struct CurveFactory<2, 2>;
template struct CurveFactory<2, 4>;