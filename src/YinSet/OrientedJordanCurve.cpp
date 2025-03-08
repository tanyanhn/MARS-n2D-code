#include "OrientedJordanCurve.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <string>

#include "Core/Curve.h"
#include "SimplicialComplex.h"
#include "YinSet/PointsLocater.h"
#include "YinSet/SegmentedRealizableSpadjor.h"
#include "YinSet/YinSet.h"

using std::string;

template <int Dim, int Order>
auto OrientedJordanCurve<Dim, Order>::locate(const Vec<Real, Dim>& p,
                                             Real tol) const -> LocateResult {
  return LocateResult(PointsLocater(tol)(
      std::vector<OrientedJordanCurve<Dim, Order>>{*this}, {p})[0]);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(const string& parameters) {
  SimplicialComplex<Vertex> kinks;
  define(parameters, kinks);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(const string& parameters,
                                             SimplicialComplex<Vertex>& kinks) {
  assert(kinks.getDimension() == -1);
  std::stringstream iss(parameters);
  define(iss, kinks);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(std::istream& iss,
                                             SimplicialComplex<Vertex>& kinks) {
  size_t nPolys;
  // iss.read((char*)&nPolys, sizeof(size_t));
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
  size_t id;
  while (iss >> id)
    kinks.insert(Simplex<Vertex>{std::initializer_list<Vertex>{{0, id}}});

  define(points, kinks);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(
    const std::vector<Vec<Real, Dim>>& points, const vector<Vertex>& indexes) {
  SimplicialComplex<Vertex> kinks;
  for (auto id : indexes) {
    kinks.insert(Simplex<Vertex>{std::initializer_list<Vertex>{id}});
  }
  define(points, kinks);
}

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(
    const std::vector<Vec<Real, Dim>>& points,
    const SimplicialComplex<Vertex>& kinks) {
  this->knots.clear();
  this->polys.clear();
  int pre = -1;
  int pos = -1;
  int last = -1;
  if (kinks.getSimplexes().empty()) {
    Curve<Dim, Order> res = fitCurve<Order>(points, Curve<2, Order>::periodic);
    this->knots = res.getKnots();
    this->polys = res.getPolys();
  } else {
    for (auto& simplex : kinks.getSimplexes()[0]) {
      assert(simplex.vertices.begin()->second < points.size() &&
             "kinks value should be points index.");
      if (pre == -1) {
        pre = simplex.vertices.begin()->second;
        last = pre;
      } else {
        pos = simplex.vertices.begin()->second;
        vector<Vec<Real, Dim>> subPoints(std::next(points.begin(), pre),
                                         std::next(points.begin(), pos + 1));
        assert(subPoints.size() > 1 && "fitCurve Point num < 2.");
        auto type = Curve<2, Order>::notAknot;
        if (subPoints.size() <= 3) type = Curve<2, Order>::nature;

        Curve<Dim, Order> res = fitCurve<Order>(subPoints, type);
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
    size_t cut;
    if (last != 0) {
      cut = this->knots.size() - 1 + subPoints.size();
      subPoints.insert(subPoints.end(), std::next(points.begin(), 1),
                       std::next(points.begin(), last + 1));
    }
    assert(subPoints.size() > 1 && "fitCurve Point num < 2.");
    Curve<Dim, Order> res = fitCurve<Order>(subPoints, Curve<2, Order>::nature);
    Real startT;
    if (this->knots.empty()) {
      startT = 0;
    } else {
      startT = this->knots.back();
      this->knots.pop_back();
      --cut;
    }
    for (auto t : res.getKnots()) {
      this->knots.push_back(startT + t);
    }
    this->polys.insert(this->polys.end(), res.getPolys().begin(),
                       res.getPolys().end());
    if (last != 0) {
      vector<Curve<Dim, Order>> curs;
      this->split(vector<Real>{this->knots[cut]}, curs, 1e-10);
      auto size = this->knots.size();
      auto startT = curs[0].getKnots()[0];
      for (size_t i = 0; i < size; ++i)
        this->knots[i] = curs[0].getKnots()[i] - startT;
      this->polys = curs[0].getPolys();
    }
  }

  // makeSelfMonotonic(distTol());
}

template <int Order>
void Circle<Order>::define(const std::string& parameters) {
  SimplicialComplex<Vertex> kinks;
  define(parameters, kinks);
}

template <int Order>
void Circle<Order>::define(const std::string& parameters,
                           SimplicialComplex<Vertex>& kinks) {
  std::stringstream iss(parameters);
  define(iss, kinks);
}

template <int Order>
void Circle<Order>::define(std::istream& iss,
                           SimplicialComplex<Vertex>& kinks) {
  assert(kinks.getDimension() == -1 && "output kinks should be empty.");
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
  for (size_t i = 0; (size_t)i <= nPolys; ++i) {
    points[i][0] = center[0] + radius * std::cos(sign * (int)i * dtheta);
    points[i][1] = center[1] + radius * std::sin(sign * (int)i * dtheta);
  }
  OrientedJordanCurve<2, Order>::define(points, kinks);
}

template <int Order>
void Rectangle<Order>::define(const std::string& parameters) {
  SimplicialComplex<Vertex> kinks;
  define(parameters, kinks);
}

template <int Order>
void Rectangle<Order>::define(const std::string& parameters,
                              SimplicialComplex<Vertex>& kinks) {
  std::stringstream iss(parameters);
  define(iss, kinks);
}

template <int Order>
void Rectangle<Order>::define(std::istream& iss,
                              SimplicialComplex<Vertex>& kinks) {
  assert(kinks.getDimension() == -1 && "output kinks should be empty.");
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
  size_t rowPolys = (bigEnd[0] - smallEnd[0]) / hL + 1;
  size_t colPolys = (bigEnd[1] - smallEnd[1]) / hL + 1;
  Real dx = (bigEnd[0] - smallEnd[0]) / rowPolys;
  Real dy = (bigEnd[1] - smallEnd[1]) / colPolys;

  std::vector<Vec<Real, 2>> points;
  vector<Vertex> kinks1;
  kinks1.push_back({0, points.size()});
  for (size_t i = 0; i < rowPolys; ++i)
    points.push_back(Vec<Real, 2>{smallEnd[0] + i * dx, smallEnd[1]});
  kinks1.push_back({0, points.size()});
  for (size_t i = 0; i < colPolys; ++i)
    points.push_back(Vec<Real, 2>{bigEnd[0], smallEnd[1] + i * dy});
  kinks1.push_back({0, points.size()});
  for (size_t i = 0; i < rowPolys; ++i)
    points.push_back(Vec<Real, 2>{bigEnd[0] - i * dx, bigEnd[1]});
  kinks1.push_back({0, points.size()});
  for (size_t i = 0; i < colPolys; ++i)
    points.push_back(Vec<Real, 2>{smallEnd[0], bigEnd[1] - i * dy});
  points.push_back(smallEnd);

  Real c = cos(theta);
  Real s = sin(theta);
  Real x;
  Real y;
  for (auto& p : points) {
    x = c * p[0] - s * p[1];
    y = s * p[0] + c * p[1];
    p[0] = x, p[1] = y;
  }
  if (!orientation) {
    std::reverse(points.begin(), points.end());
    auto size = points.size() - 1;
    for (auto& simplex : kinks1) {
      auto j = size - simplex.second;
      j = (j != size ? j : 0);
      kinks.insert(Simplex<Vertex>{std::initializer_list<Vertex>{{0, j}}});
    }
  } else {
    for (auto& vertex : kinks1) {
      kinks.insert(Simplex<Vertex>{std::initializer_list<Vertex>{vertex}});
    }
  }

  OrientedJordanCurve<2, Order>::define(points, kinks);
}

template <int Order>
std::unique_ptr<OrientedJordanCurve<2, Order>>
CurveFactory<2, Order>::createCurve(const std::string& parameters) {
  SimplicialComplex<Vertex> kinks;
  return createCurve(parameters, kinks);
}

template <int Order>
std::unique_ptr<OrientedJordanCurve<2, Order>>
CurveFactory<2, Order>::createCurve(const std::string& parameters,
                                    SimplicialComplex<Vertex>& kinks) {
  assert(kinks.getDimension() == -1);
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
  SimplicialComplex<Vertex> kinks;
  for (size_t i = 0; i < nCurves; ++i) {
    SimplicialComplex<Vertex> tmp;
    curves[i] = *createCurve(parameters[i + 1], tmp);
    if (tmp.getDimension() == 0) {
      for (auto& simplex : tmp.getSimplexes()[0]) {
        auto j = simplex.vertices.begin()->second;
        kinks.insert(Simplex<Vertex>{std::initializer_list<Vertex>{{i, j}}});
      }
    }
  }
  SegmentedRealizableSpadjor<Order> segmentedSpadjor(curves, 0);
  YinSet<2, Order> ret(segmentedSpadjor, tol);
  ret.kinks = kinks;
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