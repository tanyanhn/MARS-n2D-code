#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "Core/VecCompare.h"
#include "YinSet/OrientedJordanCurve.h"
using std::string;
using rVec = Vec<Real, 2>;

bool verifySpline(const std::vector<Polynomial<4, rVec>>& polys,
                  const std::vector<Real>& knots,
                  const bool periodic,
                  const Real tol) {
  size_t lidx, ridx;
  size_t npolys = polys.size();
  VecCompare<Real, 2> pt_cmp(tol);
  Real t = knots[npolys] - knots[npolys - 1];
  if (!periodic) {
    auto poly = polys[0];
    rVec M0 = poly.der().der()(knots[0]);
    poly = polys[npolys - 1];
    rVec Mn = poly.translate(knots[npolys - 1]).der().der()(t);
    if (pt_cmp.compare(M0, {0, 0}) != 0)
      return false;
    if (pt_cmp.compare(Mn, {0, 0}) != 0)
      return false;
  }
  for (size_t i = 0; i != polys.size(); i++) {
    if (i == 0) {
      if (!periodic)
        continue;
      lidx = polys.size() - 1;
    } else {
      lidx = i - 1;
    }
    ridx = i;
    t = knots[lidx + 1] - knots[lidx];
    // std::cout << polys[lidx](knots[lidx + 1] - knots[lidx]) <<
    // polys[lidx].translate(knots[lidx])(knots[lidx] + 1) << std::endl;
    auto lpoly = polys[lidx];
    auto rpoly = polys[ridx];
    std::cout << lpoly(t) << ' ' << rpoly[0] << std::endl;
    if (pt_cmp.compare(lpoly(t), rpoly[0]) != 0)
      ;  // return false;
    auto d_lpoly = lpoly.der();
    auto d_rpoly = rpoly.der();
    std::cout << d_lpoly(t) << d_rpoly[0] << std::endl;
    if (pt_cmp.compare(d_lpoly(t), d_rpoly[0]) != 0)
      ;  // return false;
    auto dd_lpoly = d_lpoly.der();
    auto dd_rpoly = d_rpoly.der();
    std::cout << dd_lpoly(t) << dd_rpoly[0] << std::endl;
    std::cout << std::endl;
    if (pt_cmp.compare(dd_lpoly(t), dd_rpoly[0]) != 0)
      ;  // return false;
  }
  return true;
};

template <int Order>
void drawCurve(Curve<2, Order>& cur, size_t num_piece, string of_name) {
  auto polys = cur.getPolys();
  auto knots = cur.getKnots();
  std::ofstream Curve_build_out(of_name);
  Curve_build_out << "X=[";
  Curve_build_out << polys[0][0][0] << ", " << polys[0][0][1] << ";"
                  << std::endl;
  for (size_t j = 0; j != polys.size(); j++) {
    Real leng = knots[j + 1] - knots[j];
    for (size_t i = 1; i != num_piece + 1; i++) {
      Curve_build_out << (polys[j](leng * 1.0 * i / num_piece))[0] << ", "
                      << (polys[j](leng * 1.0 * i / num_piece))[1] << ";"
                      << std::endl;
    }
  }
  Curve_build_out << "];" << std::endl;
  Curve_build_out << "hold on; plot(X(:, 1), X(:, 2), '-b');" << std::endl;
}

int main(int argc, char* argv[]) {
  string inJordanCurve1 = "OrientedJordanCurve 4 1 0 0 1 -1 0 0 -1 1 0";
  string inJordanCurve2 = "OrientedJordanCurve 4 0 0 1 0 1 2 0 2 0 0 0 1 2 3 4";
  string inJordanCurve3 =
      "OrientedJordanCurve 4 0 0 0 1 -2 1 -2 0 0 0 0 1 2 3 4";
  string inCircle1 = "Circle 0 0 1 1 1";
  string inCircle2 = "Circle 10 10 1 1 0.2";
  string inCircle3 = "Circle 15 15 1 1 0.2";
  string inRectangle1 = "Rectangle 0 0 1 2 0 1 0.5";
  string inRectangle2 =
      "Rectangle 0 0 1 2 " + std::to_string(0.5 * M_PI) + " 1 0.5";
  string inRectangle3 = "Rectangle 4 4 10 10 0 1 2";
  Real tol = 1e-8;
  CurveFactory<2, 4> factory4;
  CurveFactory<2, 2> factory2;

  auto cur = *factory4.createCurve(inCircle1);
  //std::cout << verifySpline(cur.getPolys(), cur.getKnots(), true, tol) << std::endl;
  auto polys = cur.getPolys();
  auto knots = cur.getKnots();
  drawCurve(cur, 10, "drawCircle.m");

  string input_name("../../test/data/testOrientedJordanCurve");
  int testcase;
  std::cout << "Input test case number:" << std::endl;
  std::cin >> testcase;
  // 0.one closed circle.
  // 1.unclosed circle.
  // 2.two equivalent curves.
  // 3.two intersected curves.
  std::cout << "TestCase: " << testcase << std::endl;
  std::ifstream input(input_name + std::to_string(testcase) + ".input");
  std::vector<string> factory_params(1);
  std::ostringstream oss;
  oss << std::setprecision(12) << tol;
  factory_params[0] = oss.str();
  std::string s;
  getline(input, s);
  while (!input.fail()) {
    factory_params.push_back(s);
    getline(input, s);
  }
  std::cout << "Test Order=2:" << std::endl;
  auto ys2 = factory2.createYinSet(factory_params);
  std::cout << ys2.getHasseString() << std::endl;
  std::ofstream of2(
      "result/resultOrientedJordanCurveFactory2_" + std::to_string(testcase) + ".dat",
      std::ios_base::binary);
  ys2.dump(of2);
  /*std::cout << "Test Order=4:" << std::endl;
  auto ys4 = factory4.createYinSet(factory_params);
  std::cout << ys4.getHasseString() << std::endl;
  std::ofstream of4(
      "result/resultOrientedJordanCurveFactory4_" + std::to_string(testcase) + ".dat",
      std::ios_base::binary);
  ys4.dump(of4);
  */
  // outYS(inCircle, "../result/resultCircle1.dat");
  // outYS(inRectangle2, "../result/resultRectangle2.dat");
  /*  auto compareYS = [tol, &factory2, &factory_params](string param1, string
  param2) { factory_params[1] = param1; auto ys1 =
  factory2.createYinSet(factory_params); factory_params[1] = param2; auto ys2 =
  factory2.createYinSet(factory_params);
                     // CPPUNIT_ASSERT(ys1.equal(ys2, 0.0));
                   };
  compareYS(inJordanCurve1, inCircle1);
  compareYS(inJordanCurve2, inRectangle1);
  compareYS(inJordanCurve3, inRectangle2);
  */
};
