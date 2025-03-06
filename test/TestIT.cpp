#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <vector>

#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/DIRK.h"
#include "InterfaceTracking/ERK.h"
#include "InterfaceTracking/MARS2D.h"
#include "InterfaceTracking/MARS2DIMV.h"
#include "InterfaceTracking/TestExample.h"

using namespace std;

using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;

Crv output(const Crv &crv, Point center, Real radio) {
  Vector<Point> npts;
  auto polys = crv.getPolys();
  auto knots = crv.getKnots();
  int m = 32;
  for (int i = 0; i < (int)polys.size(); i++) {
    for (int j = 0; j < m; j++) {
      Real dt = j * (knots[i + 1] - knots[i]) / m;
      auto pt = polys[i](dt);
      npts.push_back(pt);
    }
  }
  npts.push_back(polys[0][0]);
  return fitCurve<4>(npts, Curve<2, 4>::periodic);
}

void testIT() {
  int loop = 1;
  bool plot = true;
  Real tol = 1e-15;
  int stages = 5;
  cout << setiosflags(ios::scientific) << setprecision(2);
  int testCase = 10;
  TestIT test = getTest(testCase);

  // set the initial curve
  int n;
  Real dt;
  int opstride;
  Real radio = test.radio;
  Point center = test.center;
  Vector<Curve<2, 4>> crvs;
  OrientedJordanCurve<2, 4> crv;

  ERK<2, RK::ClassicRK4, VectorFunction> ERK;
  DIRK<2, RK::ESDIRK4, VectorFunction> ESDIRK4;
  DIRK<2, RK::SDIRK2, VectorFunction> SDIRK2;

  Vector<Real> time1(2 * stages - 1, 0);
  Vector<Real> time2(2 * stages - 1, 0);
  clock_t begin, end;
  string dir = std::string(ROOT_DIR) + "/results/" + test.name + "/";
  // auto success = 
  mkdir(dir.c_str(), 0755);
  // assert(success == 0);
  for (int lp = 0; lp < loop; lp++) {
    n = test.n;    //
    dt = test.dt;  //
    opstride = test.opstride;
    for (int k = 0; k < stages; k++) {
      // get the initial curve
      Vector<Point> pts;
      pts.push_back({center[0] + radio, center[1]});
      for (int i = 1; i < n; i++) {
        pts.push_back({center[0] + radio * cos(2 * M_PI / n * i),
                       center[1] + radio * sin(2 * M_PI / n * i)});
      }
      pts.push_back({center[0] + radio, center[1]});
      crv = fitCurve<4>(pts, Curve<2, 4>::periodic);
      Vector<OrientedJordanCurve<2, 4>> vcrv{crv};
      YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);

      // set the CubicMARS method
      MARS2DIMV<4, VectorFunction> CM(&ERK, 4 * M_PI * radio / n, test.rtiny);
      // MARS2DIMV<4, VectorFunction> CM(&ESDIRK4, 4 * M_PI * radio / n,
      // test.rtiny); MARS2DIMV<4, VectorFunction> CM(&SDIRK2, 4 * M_PI * radio
      // / n, test.rtiny);

      ostringstream tmps;
      tmps << k;
      string fname = dir + "No" + tmps.str();

      begin = clock();
      if (plot)
        CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname,
                          opstride);
      else
        CM.trackInterface(*test.velocity, YS, 0, dt, test.T);
      end = clock();
      time1[2 * k] += (double)(end - begin) / CLOCKS_PER_SEC;
      // get the curve after tracking
      crv = (YS.getBoundaryCycles())[0];
      crvs.push_back(crv);
      n *= 2;
      opstride *= 2;
      dt /= 2;
    }
  }

  // get the approx solution
  n *= 8;  // ensure that the chdlength is smaller than computational solutions'
  Vector<Point> rpts;
  rpts.push_back({center[0] + radio, center[1]});
  for (int i = 1; i < n; i++) {
    rpts.push_back({center[0] + radio * cos(2 * M_PI / n * i),
                    center[1] + radio * sin(2 * M_PI / n * i)});
  }
  rpts.push_back({center[0] + radio, center[1]});
  auto rcrv = fitCurve<4>(rpts, Curve<2, 4>::periodic);

  // output the convergency rate
  auto it1 = crvs.begin();
  auto it2 = crvs.begin() + stages;
  // auto it3 = crvs.end();
  auto result1 = exactError(Vector<Crv>(it1, it2), rcrv, tol);
  // auto result2 = exactError(Vector<Crv>(it1, it2), r2crv, tol);
  // auto result2 = exactError(Vector<Crv>(it2, it3), rcrv, tol);

  // auto result1 = richardsonError(crvs, tol);
  cout << "error: ";
  for (auto &i : result1) {
    cout << i << "  ";
  }
  cout << endl;
  /*
  for (auto &i : result2)
  {
      cout << i << "  ";
  }
  cout << endl;
  */

  for (auto &i : time1) {
    i = i / loop;
  }

  for (auto &i : time2) {
    i = i / loop;
  }

  for (int i = 0; i < stages - 1; i++) {
    time1[2 * i + 1] = log(time1[2 * i + 2] / time1[2 * i]) / log(2);
    time2[2 * i + 1] = log(time2[2 * i + 2] / time2[2 * i]) / log(2);
  }
  cout << "time : ";
  for (auto &i : time1) {
    cout << i << "  ";
  }
  cout << endl;
  /*
  cout << "method2 time: ";
  for (auto &i : time2)
  {
      cout << i << "  ";
  }
  cout << endl;
  */
}

void testKinks_0nk() {
  Real tol = 1e-15;
  cout << setiosflags(ios::scientific) << setprecision(2);

  // set the initial curve
  Vector<Curve<2, 4>> crvs;
  OrientedJordanCurve<2, 4> crv;

  ERK<2, RK::ClassicRK4, VectorFunction> ERK;
  DIRK<2, RK::ESDIRK4, VectorFunction> ESDIRK4;
  DIRK<2, RK::SDIRK2, VectorFunction> SDIRK2;

  // get the initial curve
  Vector<Point> pts;
  SimplicialComplex<typename YinSet<2, 4>::PointIndex> kinks;
  std::vector<typename YinSet<2, 4>::PointIndex> kps;
  for (int i = 0; i < 5; i++) {
    pts.push_back({0.2 * i, -1});
  }
  kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
      std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, 5}}});
  kps.push_back(std::make_pair(0, 5));
  for (int i = 0; i < 10; i++) {
    pts.push_back({1, -1 + 0.2 * i});
  }
  kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
      std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, 15}}});
  kps.push_back(std::make_pair(0, 15));
  for (int i = 0; i < 10; i++) {
    pts.push_back({1 - 0.2 * i, 1});
  }
  kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
      std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, 25}}});
  kps.push_back(std::make_pair(0, 25));
  for (int i = 0; i < 10; i++) {
    pts.push_back({-1, 1 - 0.2 * i});
  }
  kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
      std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, 35}}});
  kps.push_back(std::make_pair(0, 35));
  for (int i = 0; i < 5; i++) {
    pts.push_back({-1 + 0.2 * i, -1});
  }
  pts.push_back(pts[0]);
  crv.define(pts, kinks);

  Vector<OrientedJordanCurve<2, 4>> vcrv{crv};
  YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);
  YS.resetAllKinks(kps);

  // set the CubicMARS method
  MARS2DIMV<4, VectorFunction> CM(&ERK, 0.3, 0.5);
  CM.trackInterface(*squareshrink, YS, 0, 0.1, 0.4);
  // CM.trackInterface(*translation, YS, 0, 0.1, 1);
  crv = (YS.getBoundaryCycles())[0];
  auto knots = crv.getKnots();
  auto polys = crv.getPolys();
  int ln = polys.size();

  for (int i = 0; i < ln; i++) {
    cout << crv(knots[i]) << "  ";
  }
  cout << crv(knots[ln]);
  cout << endl;
}

void testKinks_0k() {
  Real tol = 1e-15;
  int stage = 4;
  std::cout << setiosflags(ios::scientific) << setprecision(2);
  ERK<2, RK::ClassicRK4, VectorFunction> ERK;
  DIRK<2, RK::ESDIRK4, VectorFunction> ESDIRK4;
  DIRK<2, RK::SDIRK2, VectorFunction> SDIRK2;
  TestIT test = getTest(9);
  int n = 80;
  Point center = test.center;
  Real h = 2.0 * test.radio / n;
  Real dt = test.dt;
  Vector<YinSet<2, 4>> vys;

  for (int k = 0; k < stage; k++) {
    // set the initial curve
    Vector<Curve<2, 4>> crvs;
    OrientedJordanCurve<2, 4> crv;
    Vector<Point> pts;
    SimplicialComplex<typename YinSet<2, 4>::PointIndex> kinks;
    std::vector<typename YinSet<2, 4>::PointIndex> kps;
    for (int i = 0; i < n; i++) {
      pts.push_back({center[0] - test.radio + h * i, center[1] - test.radio});
    }
    kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
        std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, 0}}});
    kps.push_back(std::make_pair(0, 0));
    for (int i = 0; i < n; i++) {
      pts.push_back({center[0] + test.radio, center[1] - test.radio + h * i});
    }
    kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
        std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, n}}});
    kps.push_back(std::make_pair(0, n));
    for (int i = 0; i < n; i++) {
      pts.push_back({center[0] + test.radio - h * i, center[1] + test.radio});
    }
    kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
        std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, 2 * n}}});
    kps.push_back(std::make_pair(0, 2 * n));
    for (int i = 0; i < n; i++) {
      pts.push_back({center[0] - test.radio, center[1] + test.radio - h * i});
    }
    kinks.insert(Simplex<typename YinSet<2, 4>::PointIndex>{
        std::initializer_list<typename YinSet<2, 4>::Vertex>{{0, 3 * n}}});
    kps.push_back(std::make_pair(0, 3 * n));
    pts.push_back(pts[0]);
    crv.define(pts, kinks);

    Vector<OrientedJordanCurve<2, 4>> vcrv{crv};
    YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);
    YS.resetAllKinks(kps);

    // set the CubicMARS method
    MARS2DIMV<4, VectorFunction> CM(&ESDIRK4, 2 * h, test.rtiny);
    CM.trackInterface(*test.velocity, YS, 0, dt, test.T);
    // CM.trackInterface(*translation, YS, 0, 0.1, 1);
    vys.push_back(YS);

    n *= 2;
    h /= 2;
    dt /= 2;
  }
  /*
  crv = (YS.getBoundaryCycles())[0];
  auto knots = crv.getKnots();
  auto polys = crv.getPolys();
  int ln = polys.size();

  for (int i = 0; i < ln; i++)
  {
      std::cout << crv(knots[i]) << "  ";
  }
  std::cout << crv(knots[ln]);
  std::cout << endl;
  */
  Vector<Real> result = squarerror(vys, center, test.radio);
  for (auto &i : result) {
    cout << i << "  ";
  }
  cout << endl;
}

void testKinks_circle() {
  bool plot = false;
  Real tol = 1e-15;
  int stages = 5;
  cout << setiosflags(ios::scientific) << setprecision(2);
  TestIT test = getTest(8);

  // set the initial curve
  int n;
  Real dt;
  int opstride;
  Real radio = test.radio;
  Point center = test.center;
  Vector<Curve<2, 4>> crvs;
  OrientedJordanCurve<2, 4> crv;

  ERK<2, RK::ClassicRK4, VectorFunction> ERK;
  DIRK<2, RK::ESDIRK4, VectorFunction> ESDIRK4;
  DIRK<2, RK::SDIRK2, VectorFunction> SDIRK2;

  Vector<Real> time1(2 * stages - 1, 0);
  Vector<Real> time2(2 * stages - 1, 0);
  clock_t begin, end;

  n = test.n;    //
  dt = test.dt;  //
  opstride = test.opstride;
  for (int k = 0; k < stages; k++) {
    // get the initial curve
    Vector<Point> pts;
    pts.push_back({center[0] + radio, center[1]});
    for (int i = 1; i < n; i++) {
      pts.push_back({center[0] + radio * cos(2 * M_PI / n * i),
                     center[1] + radio * sin(2 * M_PI / n * i)});
    }
    pts.push_back({center[0] + radio, center[1]});
    Vector<typename YinSet<2, 4>::PointIndex> kinks;
    kinks.push_back(std::make_pair(0, 0));
    // kinks.push_back(std::make_pair(0, n / 8));
    // kinks.push_back(std::make_pair(0, n / 4));
    // kinks.push_back(std::make_pair(0, n / 2));
    // kinks.push_back(std::make_pair(0, n * 3 / 4));
    crv = fitCurve<4>(pts, Curve<2, 4>::periodic);
    Vector<OrientedJordanCurve<2, 4>> vcrv{crv};
    YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);
    YS.resetAllKinks(kinks);

    // set the CubicMARS method
    MARS2DIMV<4, VectorFunction> CM(&SDIRK2, 4 * M_PI * radio / n, test.rtiny);

    ostringstream tmps;
    tmps << k;
    string fname = "results" + test.name + "/No" + tmps.str();

    begin = clock();

    if (plot == true)
      CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname,
                        opstride);
    else
      CM.trackInterface(*test.velocity, YS, 0, dt, test.T);

    end = clock();
    time1[2 * k] += (double)(end - begin) / CLOCKS_PER_SEC;
    // get the curve after tracking
    crv = (YS.getBoundaryCycles())[0];
    crvs.push_back(crv);
    n *= 2;
    opstride *= 2;
    dt /= 2;
  }

  // get the approx solution
  /*
  n *= 8; //ensure that the chdlength is smaller than computational solutions'
  Vector<Point> rpts;
  rpts.push_back({center[0] + radio, center[1]});
  for (int i = 1; i < n; i++)
  {
      rpts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] +
  radio * sin(2 * M_PI / n * i)});
  }
  rpts.push_back({center[0] + radio, center[1]});
  auto rcrv = fitCurve<4>(rpts, Curve<2, 4>::periodic);
  */

  // output the convergency rate
  auto result1 = circleerror(crvs, test.center, test.radio);
  for (auto &i : result1) {
    cout << i << "  ";
  }
  cout << endl;

  for (int i = 0; i < stages - 1; i++) {
    time1[2 * i + 1] = log(time1[2 * i + 2] / time1[2 * i]) / log(2);
  }
  cout << "method1 time: ";
  for (auto &i : time1) {
    cout << i << "  ";
  }
  cout << endl;
}

int main() {
  testIT();
  // testKinks_0nk();
  //  testKinks_0k();
  // testKinks_circle();
  return 0;
}
