#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#include "Core/Curve.h"
#include "Core/VecCompare.h"
#include "Core/dirConfig.h"
#include "YinSet/OrientedJordanCurve.h"
#include "YinSet/YinSet.h"
using std::string;
using rVec = Vec<Real, 2>;

bool verifySpline(const std::vector<Polynomial<4, rVec>>& polys,
                  const std::vector<Real>& knots, const bool periodic,
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
    if (pt_cmp.compare(M0, {0, 0}) != 0) return false;
    if (pt_cmp.compare(Mn, {0, 0}) != 0) return false;
  }
  for (size_t i = 0; i != polys.size(); i++) {
    if (i == 0) {
      if (!periodic) continue;
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
void drawCurve(const Curve<2, Order>& cur, size_t num_piece, string of_name) {
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
  Curve_build_out << "hold on; plot(X(:, 1), X(:, 2)); figure(1)" << std::endl;
}
template <int Order>
void drawCurve(const OrientedJordanCurve<2, Order>& cur, size_t num_piece,
               string of_name) {
  Curve<2, Order> cu = cur;
  drawCurve(cu, num_piece, of_name);
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
  // std::cout << verifySpline(cur.getPolys(), cur.getKnots(), true, tol) <<
  // std::endl;
  auto polys = cur.getPolys();
  auto knots = cur.getKnots();
  drawCurve(cur, 10, "drawCircle.m");

  string input_name(std::string(ROOT_DIR) + "/test/data/testCurveFactory");
  int testcase = 0;
  // std::cout << "Input test case number:" << std::endl;
  // std::cin >> testcase;
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
  std::ofstream of2(std::string(ROOT_DIR) +
                        "/test/result/resultOrientedJordanCurveFactory2_" +
                        std::to_string(testcase) + ".dat",
                    std::ios_base::binary);
  ys2.dump(of2);
  std::vector<std::string> params(2);
  params[0] = factory_params[0];
  params[1] = inRectangle1;
  auto ys = factory2.createYinSet(params);
  auto intersect_ys = intersect(ys, ys2, tol);
  std::cout << "Intersect test:" << std::endl;
  std::cout << intersect_ys.getHasseString() << std::endl;
  auto complement_ys = ys2.complement(tol);
  std::cout << "Complement test:" << std::endl;
  std::cout << complement_ys.getHasseString() << std::endl;
  /*std::cout << "Test Order=4:" << std::endl;
  auto ys4 = factory4.createYinSet(factory_params);
  std::cout << ys4.getHasseString() << std::endl;
  std::ofstream of4(
      "result/resultOrientedJordanCurveFactory4_" + std::to_string(testcase) +
  ".dat", std::ios_base::binary); ys4.dump(of4);
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

  std::cout << "Input kinks test case number:" << std::endl;
  // std::cin >> testcase;
  testcase = 2;
  std::cout << "kinks TestCase: " << testcase << std::endl;
  std::ifstream input2(input_name + std::to_string(testcase) + ".input");
  std::vector<string> factory_params2(1);
  factory_params2[0] = oss.str();
  getline(input2, s);
  while (!input2.fail()) {
    factory_params2.push_back(s);
    getline(input2, s);
  }
  using Vertex = YinSet<2, 4>::Vertex;
  std::cout << "Test Order=4:" << std::endl;
  auto ys4 = factory4.createYinSet(factory_params2);
  auto sims = ys4.getKinks();
  assert(sims.getSimplexes()[0].size() == 8);
  assert(sims.contain(Simplex<Vertex>{{{1ul, 3ul}}}));
  drawCurve(ys4.getBoundaryCycles()[0], 10, "kinks0.m");
  ys4.resetAllKinks(vector<typename YinSet<2, 4>::PointIndex>());
  sims = ys4.getKinks();
  assert(sims.getSimplexes().size() == 0);
  drawCurve(ys4.getBoundaryCycles()[0], 10, "kinks1.m");
  ys4.insertKink(std::make_pair(0, 8));
  sims = ys4.getKinks();
  assert(sims.getSimplexes()[0].size() == 1);
  assert(sims.contain(Simplex<Vertex>{{{0ul, 8ul}}}));
  drawCurve(ys4.getBoundaryCycles()[0], 10, "kinks2.m");
  ys4.eraseKink({0, 8});
  sims = ys4.getKinks();
  assert(sims.getSimplexes().size() == 0);
  drawCurve(ys4.getBoundaryCycles()[0], 10, "kinks3.m");
  ys4.resetAllKinks(
      vector<typename YinSet<2, 4>::PointIndex>({{0, 0}, {0, 8}, {1, 9}}));
  sims = ys4.getKinks();
  assert(sims.getSimplexes()[0].size() == 3);
  assert(sims.contain(Simplex<Vertex>{{{1ul, 9ul}}}));
  drawCurve(ys4.getBoundaryCycles()[0], 10, "kinks40.m");
  drawCurve(ys4.getBoundaryCycles()[1], 10, "kinks41.m");
  drawCurve(ys4.getBoundaryCycles()[2], 10, "kinks50.m");
  ys4.resetAllKinks(
      vector<typename YinSet<2, 4>::PointIndex>({{2, 2}, {2, 3}, {2, 4}}));
  drawCurve(ys4.getBoundaryCycles()[2], 10, "kinks51.m");
  ys4.insertKink(std::make_pair(2, 0));
  drawCurve(ys4.getBoundaryCycles()[2], 10, "kinks52.m");
};
