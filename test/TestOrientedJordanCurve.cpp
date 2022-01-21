#include "TestOrientedJordanCurve.h"
#include <string>
#include "YinSet/OrientedJordanCurve.h"

using std::string;

void TestOrientedJordanCurve::testOrientedJordanCurve(const std::string& input, const Real tol)
{
  VecCompare<Real, 2> pt_cmp(tol);
  std::stringstream iss(parameters);
  unsigned int nPolys;
  iss >> nPolys;
  std::vector<rVec> points(nPolys + 1);
  for(size_t i = 0; i != points.size(); i++) {
    iss >> points[i][0] >> points[i][1];
  }
  std::vector<unsigned int> kink_ids;
  unsigned int idx;
  while(iss >> id)
    kinkid.push_back(id);
  //case Order = 2
  {
    OrientedJordanCurve<2, 2> jordancurve2;
    jordancurve2.define(input);
    auto knots = jordancurve2.getKnots();
    auto polys = jordancurve2.getPolys();
    Real t = knots.back() - knots[knots.size() - 2];
    //verify knots
    for(size_t i = 0; i != polys.size(); i++) {
      CPPUNIT_ASSERT(pt_cmp.compare(points[i], polys[i][0]) == 0);
    }
    CPPUNIT_ASSERT(pt_cmp.compare(points.back(), polys.back()(t)) == 0);
  }
  //case Order = 4
  {
    OrientedJordanCurve<2, 4> jordancurve4;
    jordancurve4.define(input);
    auto knots = jordancurve4.getKnots();
    auto polys = jordancurve4.getPolys();
    Real t = knots.back() - knots[knots.size() - 2];
    //verify knots
    for(size_t i = 0; i != polys.size(); i++) {
      CPPUNIT_ASSERT(pt_cmp.compare(points[i], polys[i][0]) == 0);
    }
    CPPUNIT_ASSERT(pt_cmp.compare(points.back(), polys.back()(t)) == 0);
    //verify polys
    if(kink_ids.size() == 0)
      CPPUNIT_ASSERT(verifySpline(polys, knots, true, tol));
    else {
      int pre = -1, pos = -1;
      for(idx : kink_ids) {
        if(pre == -1) {
          pre = idx;
          continue;
        }
        pos = idx;
        vector<Polynomial<4, rVec>> subpolys(std::next(polys.begin(), pre),
                                             std::next(polys.begin(), pos));
        vector<Real> subknots(std::next(knots.begin(), pre),
                              std::next(knots.begin(), pos + 1));
        CPPUNIT_ASSERT(verifySpline(subpolys, knots, false, tol));
        pre = pos;
      }
    }
  }
};

void TestOrientedJordanCurve::testCircle(const string& input, const Real tol)
{
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
  //verify orientation
  rVec dir = polys[0](knots[1]) - polys[0][0];
  CPPUNIT_ASSERT((dir[1] > 0) == (orientation));
  //verify knots
  for(size_t i = 0; i != polys.size(); i++) {
    rVec knot2center = polys[i][0] - center;
    CPPUNIT_ASSERT(std::abs(norm(knot2center) - radius) < tol);
  }
  //verify polynomials
  CPPUNIT_ASSERT(verifySpline(polys, knots, true, tol));
};

void TestOrientedJordanCurve::testRectangle(const string& input, const Real tol)
{
  std::stringstream iss(parameters);
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
  //verify orientation
  rVec dir = polys[0](knots[1]) - polys[0][0];
  dir = normalize(dir);
  if(orientation)
    CPPUNIT_ASSERT(pt_cmp.compare(dir, {c, s}) == 0);
  else
    CPPUNIT_ASSERT(pt_cmp.compare(dir, {-s, c}) == 0);
  //verify knots
  rVec sw{smallEnd}, se{bigEnd[0],smallEnd[1]};
  rVec nw{smallEnd[0], bigEnd[1]}, ne{bigEnd};
  std::vector<rVec> corners{sw, se, nw, ne, sw};
  int pole = 0;
  for(size_t i = 0; i != polys.size(); i++) {
    Real t = knots[i + 1] - knots[i];
    rVec dir = poly[i](t) - poly[i][0];
    dir = normalize(dir);
    if(pt_cmp(poly[i][0], corners[pole]) == 0)
      pole++;
    CPPUNIT_ASSERT(pole < 5);
    auto isOrthogonal = [tol](rVec v1, rVec v2) {
                          real res = v1[0]*v2[1] + v1[1]*v2[0];
                          if(std::abs(res) < tol)
                            return true;
                          return false;
                        };
    if((pole + 1)%2 == 0) {
      CPPUNIT_ASSERT(isOrthogonal(dir, {-s, c}));
    }
    else {
      CPPUNIT_ASSERT(isOrthogonal(dir, {c, s}));
    }
  }
  //verify polys
  for(size_t i = 0; i != polys.size(); i++) {
    auto dd_poly = poly[i].der().der();
    CPPUNIT_ASSERT(pt_cmp(dd_poly[0], {0.0, 0.0}) == 0);
  }
};

bool TestOrientedJordanCurve::verifySpline(const std::vector<Polynomial<4, rVec>>& polys,
                                           const std::vector<Real>& knots,
                                           const bool periodic,
                                           const Real tol)
{
  size_t lidx, ridx;
  size_t npolys = polys.size();
  VecCompare<Real, 2> pt_cmp(tol);
  Real t = knots[npolys] - knots[npolys - 1];
  if(!periodic) {
    rVec M0 = polys[0].der().der()[0];
    rVec Mn = polys[npolys - 1].der().der()(t);
    if(pt_cmp.compare(M0, {0, 0}) != 0)
      return false;
    if(pt_cmp.compare(Mn, {0, 0}) != 0)
      return false;
  }
  for(size_t i = 0; i != polys.size(); i++) {
    if(i == 0) {
      if(!periodic)
        continue;
      lidx = polys.size() - 1;
    }
    else {
      lidx = i - 1;
    }
    ridx = i;
    Polynomial<4, rVec> lpoly, rpoly;
    t = knots[lidx + 1] - knots[lidx];
    lpoly = polys[lidx];
    rpoly = polys[ridx];
    if(pt_cmp.compare(lpoly(t), rpoly[0]) != 0)
      return false;
    auto d_lpoly = lpoly.der();
    auto d_rpoly = rpoly.der();
    if(pt_cmp.compare(d_lpoly(t), d_rpoly[0]) != 0)
      return false;
    auto dd_lpoly = d_lpoly.der();
    auto dd_rpoly = d_rpoly.der();
    if(pt_cmp.compare(dd_lpoly(t), dd_rpoly[0]) != 0)
      return false;
  }
  return true;
};

void TestOrientedJordanCurve::doTest() {
  string inJordanCurve1 = "4 1 0 0 1 -1 0 0 -1 1 0";
  string inJordanCurve2 = "4 0 0 1 0 1 2 0 2 0 0 0 1 2 3 4";
  string inJordanCurve3 = "4 0 0 0 1 -2 1 -2 0 0 0 0 1 2 3 4";
  string inCircle1 = "0 0 1 1 2";
  string inRectangle1 = "0 0 1 2 0 1 0.5";
  string inRectangle2 = "0 0 1 2 " + std::to_string(0.5*M_PI) + " 1 0.5";
  Real tol = 1e-6;
  std::vector<string> factory_params{std::to_string{tol}, {}};
  testOrientedJordanCurve(inJordanCurve1, tol);
  testCircle(inCircle1, tol);
  testRectangle(inRectangle1, tol);

  CurveFactory<2, 4> factory4;
  factory_params[2] = inJordanCurve1;
  YinSet<2, 4> YS1, YS2;
  auto compareYS = [&factory4, &factory_params](string param1, string param2) {
                     factory_params[1] = param1;
                     auto ys1 = factory4(factory_params);
                     factory_params[1] = param2;
                     auto ys2 = factory4(factory_params);
                     CPPUNIT_ASSERT(ys1.equal(ys2, tol));
                   };
  compareYS(inJordanCurve1, inCircle1);
  compareYS(inJordanCurve2, inRectangle1);
  compareYS(inJordanCurve3, inRectangle2);
};

int main() {
  TestOrientedJordanCurve test1;
  test1.doTest();
}
