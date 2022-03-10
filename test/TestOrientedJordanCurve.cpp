#include "TestOrientedJordanCurve.h"
#include <fstream>
#include <string>
#include "Core/VecCompare.h"

using std::string;

bool TestOrientedJordanCurve::verifySpline(
    const std::vector<Polynomial<4, rVec>>& polys,
    const std::vector<Real>& knots,
    const bool periodic,
    const Real tol) {
  size_t lidx, ridx;
  size_t npolys = polys.size();
  VecCompare<Real, 2> pt_cmp(tol);
  Real t = knots[npolys] - knots[npolys - 1];
  if (!periodic) {
    rVec M0 = polys[0].der().der()[0];
    rVec Mn = polys[npolys - 1].der().der()(t);
    if (pt_cmp.compare(M0, {0, 0}) != 0)
      return false;
    if (pt_cmp.compare(Mn, {0, 0}) != 0)
      return false;
  }
  //!!! unverify the start knot.
  for (size_t i = 0; i != polys.size(); i++) {
    if (i == 0) {
      if (!periodic)
        continue;
      lidx = polys.size() - 1;
    } else {
      lidx = i - 1;
    }
    ridx = i;
    Polynomial<4, rVec> lpoly, rpoly;
    t = knots[lidx + 1] - knots[lidx];
    lpoly = polys[lidx];
    rpoly = polys[ridx];
    if (pt_cmp.compare(lpoly(t), rpoly[0]) != 0)
      return false;
    auto d_lpoly = lpoly.der();
    auto d_rpoly = rpoly.der();
    if (pt_cmp.compare(d_lpoly(t), d_rpoly[0]) != 0)
      return false;
    auto dd_lpoly = d_lpoly.der();
    auto dd_rpoly = d_rpoly.der();
    if (pt_cmp.compare(dd_lpoly(t), dd_rpoly[0]) != 0)
      return false;
  }
  return true;
};

bool TestOrientedJordanCurve::testOrientedJordanCurve(const string& input,
                                                      const Real tol,
                                                      string& message) {
  VecCompare<Real, 2> pt_cmp(tol);
  std::stringstream iss(input);
  unsigned int nPolys;
  iss >> nPolys;
  std::vector<rVec> init_points(nPolys + 1);
  for (size_t i = 0; i != init_points.size(); i++) {
    iss >> init_points[i][0] >> init_points[i][1];
  }
  std::vector<unsigned int> kink_ids;
  unsigned int idx;
  while (iss >> idx)
    kink_ids.push_back(idx);
  std::vector<rVec> points(nPolys + 1);
  for (size_t i = 0; i != points.size(); i++) {
    points[i] = init_points[i];
  }
  // case Order = 2
  {
    OrientedJordanCurve<2, 2> jordancurve2;
    jordancurve2.define(input);
    auto knots = jordancurve2.getKnots();
    auto polys = jordancurve2.getPolys();
    Real t = knots.back() - knots[knots.size() - 2];
    // verify knots
    for (size_t i = 0; i != polys.size(); i++) {
      if (pt_cmp.compare(points[i], polys[i][0]) != 0) {
        message = "Load wrong knot.";
        return false;
      }
    }
    if (pt_cmp.compare(points.back(), polys.back()(t)) != 0) {
      message = "Load wrong knot.";
      return false;
    }
  }
  // case Order = 4
  {
    OrientedJordanCurve<2, 4> jordancurve4;
    jordancurve4.define(input);
    auto knots = jordancurve4.getKnots();
    auto polys = jordancurve4.getPolys();
    Real t = knots.back() - knots[knots.size() - 2];
    // verify knots
    for (size_t i = 0; i != polys.size(); i++) {
      if (pt_cmp.compare(points[i], polys[i][0]) != 0) {
        message = "Load wrong knot.";
        return false;
      }
    }
    if (pt_cmp.compare(points.back(), polys.back()(t)) != 0) {
      message = "Load wrong knot.";
      return false;
    }
    // verify polys
    if (kink_ids.size() == 0) {
      if (!verifySpline(polys, knots, true, tol)) {
        message = "Wrong spline.";
        return false;
      }
    } else {
      kink_ids.push_back(points.size() - 1);
      int pre = -1, pos = -1, first = kink_ids[0];
      for (auto index : kink_ids) {
        if (pre == -1) {
          pre = index;
          continue;
        }
        pos = index;
        vector<Polynomial<4, rVec>> subpolys(std::next(polys.begin(), pre),
                                             std::next(polys.begin(), pos));
        vector<Real> subknots(std::next(knots.begin(), pre),
                              std::next(knots.begin(), pos + 1));
        if (pos == (int)points.size() - 1) {
          subpolys.insert(subpolys.end(), std::next(polys.begin(), 0),
                          std::next(polys.begin(), first));
          Real t = subknots.back();
          for (auto it = std::next(knots.begin(), 1);
               it != std::next(knots.begin(), first + 1); ++it)
            subknots.push_back(*it + t);
        }
        if (!verifySpline(subpolys, subknots, false, tol)) {
          message = "Wrong spline.";
          return false;
        }
        pre = pos;
      }
    }
  }
  return true;
};

bool TestOrientedJordanCurve::testCircle(const string& input,
                                         const Real tol,
                                         string& message) {
  std::stringstream iss(input);
  rVec center;
  Real radius;
  bool orientation;
  Real hL;
  iss >> center[0] >> center[1] >> radius >> orientation >> hL;
  Circle<4> circle4;
  circle4.define(input);
  auto knots = circle4.getKnots();
  auto polys = circle4.getPolys();
  // verify orientation
  rVec dir = polys[0](knots[1]) - polys[0][0];
  if ((dir[1] > 0) != (orientation)) {
    message = "Wrong orientation";
    return false;
  }
  // verify knots
  for (size_t i = 0; i != polys.size(); i++) {
    rVec knot2center = polys[i][0] - center;
    if (std::abs(norm(knot2center) - radius) > tol) {
      message = "Load wrong knot.";
      return false;
    }
  }
  // verify polynomials
  if (!verifySpline(polys, knots, true, tol)) {
    message = "Wrong spline.";
    return false;
  }
  return true;
};

bool TestOrientedJordanCurve::testRectangle(const string& input,
                                            const Real tol,
                                            string& message) {
  std::stringstream iss(input);
  VecCompare<Real, 2> pt_cmp(tol);
  Vec<Real, 2> smallEnd, bigEnd;
  Real theta;
  bool orientation;
  Real hL;
  iss >> smallEnd[0] >> smallEnd[1] >> bigEnd[0] >> bigEnd[1] >> theta >>
      orientation >> hL;
  Real c = cos(theta), s = sin(theta);
  Rectangle<4> rect4;
  rect4.define(input);
  auto knots = rect4.getKnots();
  auto polys = rect4.getPolys();
  // verify orientation
  rVec dir = polys[0](knots[1]) - polys[0][0];
  dir = normalize(dir);
  if ((orientation && pt_cmp.compare(dir, {c, s}) != 0) ||
      (!orientation && pt_cmp.compare(dir, {-s, c}) != 0)) {
    message = "Wrong orientation";
    return false;
  };
  // verify knots
  rVec sw{smallEnd}, se{bigEnd[0], smallEnd[1]};
  rVec nw{smallEnd[0], bigEnd[1]}, ne{bigEnd};
  std::vector<rVec> corners{sw, se, ne, nw, sw};
  if (!orientation) {
    corners[1] = nw;
    corners[3] = se;
  }
  for (size_t i = 0; i != 5; i++) {
    Real x = c * corners[i][0] - s * corners[i][1];
    Real y = s * corners[i][0] + c * corners[i][1];
    corners[i][0] = x;
    corners[i][1] = y;
  }
  int pole = 0;
  rVec end = polys[polys.size() - 1](knots.back() - knots[knots.size() - 2]);
  if (pt_cmp.compare(end, sw) != 0) {
    message = "Load wrong knot.";
    return false;
  }
  auto isOrthogonal = [tol](rVec v1, rVec v2) {
    Real res = v1[0] * v2[0] + v1[1] * v2[1];
    if (std::abs(res) < tol)
      return true;
    return false;
  };
  for (size_t i = 0; i != polys.size(); i++) {
    Real t = knots[i + 1] - knots[i];
    rVec dir = polys[i](t) - polys[i][0];
    dir = normalize(dir);
    if (pt_cmp.compare(polys[i][0], corners[pole]) == 0)
      pole++;
    if (pole >= 5) {
      message = "Reached too many corners.";
      return false;
    }
    if ((pole + 1) % 2 == 0) {
      rVec e = orientation ? rVec{-s, c} : rVec{c, s};
      if (!isOrthogonal(dir, e)) {
        message = "Load wrong knot.";
        return false;
      }
    } else {
      rVec e = orientation ? rVec{c, s} : rVec{-s, c};
      if (!isOrthogonal(dir, e)) {
        message = "Load wrong knot.";
        return false;
      }
    }
  }
  if (pole < 4) {
    message = "Reached too few corners.";
    return false;
  }
  // verify polys
  for (size_t i = 0; i != polys.size(); i++) {
    auto dd_poly = polys[i].der().der();
    if (pt_cmp(dd_poly[0], {0.0, 0.0}) != 0) {
      message = "Wrong spline.";
      return false;
    }
  }
  return true;
};

void TestOrientedJordanCurve::doTest(int num) {
  string input_name("data/testOrientedJordanCurve-" + std::to_string(num) +
                    ".input");
  std::ifstream input(input_name);
  std::vector<string> testcase(3);
  for (int i = 0; i != 3; i++) {
    getline(input, testcase[i]);
  }
  Real tol = 1e-6;
  string message;
  CPPUNIT_ASSERT_MESSAGE(message,
                         testOrientedJordanCurve(testcase[0], tol, message));
  CPPUNIT_ASSERT_MESSAGE(message, testCircle(testcase[1], tol, message));
  CPPUNIT_ASSERT_MESSAGE(message, testRectangle(testcase[2], tol, message));
};

/*
int main() {
  TestOrientedJordanCurve test1;
  test1.doTest();
}

*/
