#include <iostream>
#include <fstream>
#include "YinSet/OrientedJordanCurve.h"
#include "Core/VecCompare.h"
using std::string;
using rVec = Vec<Real, 2>;

bool verifySpline(const std::vector<Polynomial<4, rVec>>& polys,
                  const std::vector<Real>& knots,
                  const bool periodic,
                  const Real tol)
{
  size_t lidx, ridx;
  size_t npolys = polys.size();
  VecCompare<Real, 2> pt_cmp(tol);
  Real t = knots[npolys] - knots[npolys - 1];
  if(!periodic) {
    auto poly = polys[0];
    rVec M0 = poly.der().der()(knots[0]);
    poly = polys[npolys - 1];
    rVec Mn = poly.translate(knots[npolys - 1]).der().der()(t);
    if(pt_cmp.compare(M0, {0, 0}) != 0)
      return false;
    if(pt_cmp.compare(Mn, {0, 0}) != 0)
      return false;
  }
  for(size_t i = 1; i != polys.size(); i++) {
    if(i == 0) {
      if(!periodic)
        continue;
      lidx = polys.size() - 1;
    }
    else {
      lidx = i - 1;
    }
    ridx = i;
    t = knots[lidx + 1] - knots[lidx];
    // std::cout << polys[lidx](knots[lidx + 1] - knots[lidx]) << polys[lidx].translate(knots[lidx])(knots[lidx] + 1) << std::endl;
    auto lpoly = polys[lidx];
    auto rpoly = polys[ridx];
    std::cout << lpoly(t) << ' ' << rpoly[0] << std::endl;
    if(pt_cmp.compare(lpoly(t), rpoly[0]) != 0)
      return false;
    auto d_lpoly = lpoly.der();
    auto d_rpoly = rpoly.der();
    std::cout << d_lpoly(t) << d_rpoly[0] << std::endl;
    if(pt_cmp.compare(d_lpoly(t), d_rpoly[0]) != 0)
      return false;
    auto dd_lpoly = d_lpoly.der();
    auto dd_rpoly = d_rpoly.der();
    std::cout << dd_lpoly(t) << dd_rpoly[0] << std::endl;
    std::cout << std::endl;
    if(pt_cmp.compare(dd_lpoly(t), dd_rpoly[0]) != 0)
      return false;
  }
  return true;
};

int main(int argc, char* argv[])
{
  string inJordanCurve1 = "4 1 0 0 1 -1 0 0 -1 1 0";
  string inJordanCurve2 = "4 0 0 1 0 1 2 0 2 0 0 0 1 2 3 4";
  string inJordanCurve3 = "4 0 0 0 1 -2 1 -2 0 0 0 0 1 2 3 4";
  string inCircle1 = "0 0 1 1 2";
  string inRectangle1 = "0 0 1 2 0 1 0.5";
  string inRectangle2 = "0 0 1 2 " + std::to_string(0.5*M_PI) + " 1 0.5";
  Real tol = 1e-8;
  std::ofstream Curve_rebuild_out("output_test.m");
  std::vector<string> factory_params(2);
  factory_params[0] = std::to_string(0.0);
  //  testOrientedJordanCurve(inJordanCurve1, tol);
  //testCircle(inCircle1, tol);
  //testRectangle(inRectangle1, tol);

  CurveFactory<2, 4> factory4;
  CurveFactory<2, 2> factory2;
  inJordanCurve1 = "OrientedJordanCurve " + inJordanCurve1;
  inJordanCurve2 = "OrientedJordanCurve " + inJordanCurve2;
  inJordanCurve3 = "OrientedJordanCurve " + inJordanCurve3;
  inCircle1 = "Circle " + inCircle1;
  inRectangle1 = "Rectangle " + inRectangle1;
  inRectangle2 = "Rectangle " + inRectangle2;
  auto cur = *factory4.createCurve(inCircle1);
  std::cout << verifySpline(cur.getPolys(), cur.getKnots(), true, tol) << std::endl;
  auto polys = cur.getPolys();
  auto knots = cur.getKnots();
  /*Curve_rebuild_out << "X=[";
  Curve_rebuild_out << polys[0][0][0] << ", " << polys[0][0][1] << ";" << std::endl;
  for(size_t j = 0; j != polys.size(); j++) {
    Curve_rebuild_out << (polys[j](knots[j + 1] - knots[j]))[0] << ", " << (polys[j](knots[j + 1] - knots[j]))[1] << ";" << std::endl;
  }
  Curve_rebuild_out << "];" << std::endl;
  Curve_rebuild_out << "hold on; plot(X(:, 1), X(:, 2), '-b');" << std::endl;*/

  //factory_params.push_back(inRectangle2);
  auto outYS = [&factory4, &factory_params](string param, string of_name) {
                 Curve<2, 4> curv = *factory4.createCurve(param);
                 auto knots = curv.getKnots();
                 auto polys = curv.getPolys();
                 for(auto p : polys) {
                   std::cout << p[0] << std::endl;
                 }
                 std::cout << polys.back()(knots.back() - knots[knots.size() - 2]) << '\n' << '\n';
                 factory_params[1] = param;
                 auto ys = factory4.createYinSet(factory_params);
                 std::cout << ys.getHasseString() << std::endl;
                 std::ofstream of(of_name, std::ios_base::binary);
                 ys.dump(of);
               };
  outYS(inJordanCurve1, "../result/resultOrientedJordanCurve1.dat");
  //outYS(inCircle, "../result/resultCircle1.dat");
  //outYS(inRectangle2, "../result/resultRectangle2.dat");
  auto compareYS = [tol, &factory2, &factory_params](string param1, string param2) {
                     factory_params[1] = param1;
                     auto ys1 = factory2.createYinSet(factory_params);
                     factory_params[1] = param2;
                     auto ys2 = factory2.createYinSet(factory_params);
                     // CPPUNIT_ASSERT(ys1.equal(ys2, 0.0));
                   };
  compareYS(inJordanCurve1, inCircle1);
  compareYS(inJordanCurve2, inRectangle1);
  compareYS(inJordanCurve3, inRectangle2);
};
