#include "OrientedJordanCurve.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <sstream>
#include <string>
#include "Core/Curve.h"
#include "SimplicialComplex.h"
#include "YinSet.h"
#include "YinSet/SegmentedRealizableSpadjor.h"

using std::string;

template <int Dim, int Order>
void OrientedJordanCurve<Dim, Order>::define(const string& parameters) {
  std::stringstream iss(parameters);
  unsigned int nPolys;
  // iss.read((char*)&nPolys, sizeof(unsigned int));
  iss >> nPolys;
  std::vector<Vec<Real, Dim>> points(nPolys + 1);
  // std::vector<unsigned int> kinks;
  SimplicialComplex kinks;
  for (int i = 0; i <= nPolys; i++) {
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
    Curve<Dim, Order> res = fitCurve<Order>(points);
    this->knots = res.getKnots();
    this->polys = res.getPolys();
  } else {
    int pre = -1, pos = -1;
    for (auto& simplex : kinks.getSimplexes()[0]) {
      if (pre == -1) {
        pre = *simplex.vertices.begin();
      } else {
        pos = *simplex.vertices.begin();
        vector<Vec<Real, Dim>> subPoints(std::next(points.begin(), pre),
                                         std::next(points.begin(), pos + 1));
        // TODO need use new fitCurve for natural boundary condition.
        Curve<Dim, Order> res = fitCurve<Order>(subPoints);
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
  }
}

template <int Order>
void Circle<Order>::define(const std::string& parameters) {
  std::stringstream iss(parameters);
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
  int nPolys = radius * 2 * M_PI / hL + 1;
  Real dtheta = 2 * M_PI / nPolys;
  std::vector<Vec<Real, 2>> points(nPolys + 1);
  for (size_t i = 0; i <= nPolys; ++i) {
    points[i][0] = center[0] + radius * std::cos(sign * i * dtheta);
    points[i][1] = center[1] + radius * std::sin(sign * i * dtheta);
  }
  OrientedJordanCurve<2, Order>::define(points, SimplicialComplex());
}

template <int Order>
void Rectangle<Order>::define(const std::string& parameters) {
  std::stringstream iss(parameters);
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
  int rowPolys = (bigEnd[0] - smallEnd[0]) / hL + 1,
      colPolys = (bigEnd[1] - smallEnd[1]) / hL + 1;
  Real dx = (bigEnd[0] - smallEnd[0]) / rowPolys,
       dy = (bigEnd[1] - smallEnd[1]) / colPolys;

  std::vector<Vec<Real, 2>> points;
  SimplicialComplex kinks;
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
  std::stringstream iss(parameters);
  string type;
  iss >> type;
  std::unique_ptr<OrientedJordanCurve<2, Order>> ret(nullptr);
  if (type == "OrientedJordanCurve") {
    ret.reset(new OrientedJordanCurve<2, Order>());
    ret->define(parameters.substr(20));
  } else if (type == "Circle") {
    ret.reset(new Circle<Order>());
    ret->define(parameters.substr(7));
  } else if (type == "Rectangle") {
    ret.reset(new Rectangle<Order>());
    ret->define(parameters.substr(10));
  } else {
    assert(false && "undefined type.");
  }
  return ret;
}

template <int Order>
YinSet<2, Order> CurveFactory<2, Order>::createYinSet(
    const std::vector<std::string>& parameters) {
  size_t nCurves = parameters.size() - 1;
  std::vector<Curve<2, Order>> curves(nCurves);
  Real tol = stod(parameters[0]);
  for (size_t i = 0; i < nCurves; ++i) {
    curves[i] = *createCurve(parameters[i + 1]);
  }
  SegmentedRealizableSpadjor<Order> segmentedSpadjor(curves, tol);

  return YinSet<2, Order>(segmentedSpadjor, tol);
}

template class OrientedJordanCurve<2, 2>;
template class OrientedJordanCurve<2, 4>;
template class Circle<2>;
template class Circle<4>;
template class Rectangle<2>;
template class Rectangle<4>;
template class CurveFactory<2, 2>;
template class CurveFactory<2, 4>;