#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include "GeometricFlow/DiscreteVecCurvature.h"
#include "GeometricFlow/DiscreteVecCurvatureLaplacian.h"
#include "InterfaceTracking/MARS2D.h"
#include "InterfaceTracking/MARS2DIMV.h"
#include "InterfaceTracking/ERK.h"
#include "InterfaceTracking/DIRK.h"
#include "InterfaceTracking/TestExample.h"

using namespace std;

using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;



void testFL()
{
    cout << setiosflags(ios::scientific) << setprecision(4);

    TestIT test = getTest(6);
    //set the initial curve
    int n;
    int loop = 3;
    Real dt;
    dt = test.dt;
    Real tol = 1e-15;
    Real radio = test.radio;
    Point center = test.center;
    Vector<Curve<2, 4>> crvs;
    Curve<2, 4> crv;

    DIRK<2, RK::ESDIRK4> ESDIRK4;
    ERK<2, RK::ClassicRK4> ERK;
    
    //clock_t begin, end;
    radio = 1;
    center = Point{0.0,0.0};
    n = 32;
    dt = 0.005;
    cout << "Circle test:" << endl;
    for (int j = 0 ; j < loop ; j++){
      //get the initial curve
      Vector<Point> pts;
      pts.push_back({center[0] + radio, center[1]});
      for (int i = 1; i < n; i++)
        {
          pts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
        }
      pts.push_back({center[0] + radio, center[1]});
      cout << "----------- n = " << n << ": ------------" << endl;
      crv = fitCurve<4>(pts, true);
      Vector<Curve<2, 4>> vcrv{crv};
      YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);

      //MARS2DIMV<4> CM(&ESDIRK4, 4 * M_PI * radio / n, test.rtiny);
      
      DiscreteVecCurvature<2,4> SFV1;
      DiscreteVecCurvatureLaplacian<2,2> SFV2;
      
      Vector<Point> velocity = SFV1(pts,0);
      for (auto it = velocity.begin(); it != velocity.end() ; ++it)
        cout << *it << endl;
      // Tensor<Real,2> Jacobi = SFV2.getJacobi(pts,0);
      // for (int l = 0; l < (Jacobi.size())[0]; l++){
      //   for (int m = 0; m < (Jacobi.size())[1]; m++)
      //     cout << Jacobi(l,m) << " ";
      //   cout << endl;
      // }
      //CM.trackInterface(SFV1, YS, 0, dt, 0.5);
      n*=2;
      dt /= 2;
    }
    
    n = 128;
    dt = 0.05;
    radio = 0.3;
    center = Point{1.0,1.0};
    cout << "Star test:" << endl;
    for (int j = 0 ; j < loop ; j++){
      //get the initial curve
      Vector<Point> pts;
      pts.push_back({center[0] + radio, 0});
      for (int i = 1; i < n; i++)
        {
          pts.push_back({(center[0] + radio * cos(6 * 2 * M_PI / n * i))*cos(2 * M_PI / n * i),
                         (center[1] + radio * cos(6 * 2 * M_PI / n * i))*sin(2 * M_PI / n * i)});
        }
      pts.push_back({center[0] + radio, 0});
      cout << "----------- n = " << n << ": ------------" << endl;
      crv = fitCurve<4>(pts, true);
      Vector<Curve<2, 4>> vcrv{crv};
      YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);

      //MARS2DIMV<4> CM(&ERK, 4 * M_PI * radio / n, test.rtiny);
      
      
      
      
      // DiscreteVecCurvature<2,2> SFV1;
      DiscreteVecCurvatureLaplacian<2,4> SFV2;
      
      Vector<Point> velocity = SFV2(pts,0);
      for (auto it = velocity.begin(); it != velocity.end() ; ++it)
        cout << *it << endl;
      // Tensor<Real,2> Jacobi = SFV2.getJacobi(pts,0);
      // for (int l = 0; l < (Jacobi.size())[0]; l++){
      //   for (int m = 0; m < (Jacobi.size())[1]; m++)
      //     cout << Jacobi(l,m) << " ";
      //   cout << endl;
      // }
      //CM.trackInterface(SFV1, YS, 0, dt, test.T);
      n*=2;
      dt /= 2;
    
    }


   
      //Vector<Real> arclength2 = calArcLength<3>(crv);
      // Tensor<Real,2> FD2coes2 = calFD2coes<4>(arclength2);
      // cout << "cal FD2 coes with order 4 :" << endl;
      // for (int k = 0 ; k < (int)arclength2.size()-1; k++){
      // for (int l = 0 ; l < 6 ; l++)
      //   cout << FD2coes2(l,k) << " ";
      // cout << endl;
      // }







}

int main()
{
    testFL();

    return 0;
}
