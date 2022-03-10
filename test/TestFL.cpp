#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include "GeometricFlow/VectorForCurvatureFlow.h"
#include "GeometricFlow/VectorForSurfaceDiffusionFlow.h"
#include "GeometricFlow/VectorForSurfaceDiffusionFlowNew.h"
#include "GeometricFlow/VelocityTest.h"
#include "InterfaceTracking/MARS2D.h"
#include "InterfaceTracking/MARS2DIMV.h"
#include "InterfaceTracking/ERK.h"
#include "InterfaceTracking/DIRK.h"
#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/TestExample.h"

using namespace std;

using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;



void testFL()
{
    cout << setiosflags(ios::scientific) << setprecision(4);

    TestIT test = getTest(6);
    int n;
    int loop = 1;
    Real dt;
    dt = test.dt;
    Real tol = 1e-15;
    Real radio = test.radio;
    Point center = test.center;
    Vector<Curve<2, 4>> crvs;
    Curve<2, 4> crv;

    DIRK<2, RK::ESDIRK4, VectorOnHypersurface> ESDIRK4;
    ERK<2, RK::ClassicRK4, VectorOnHypersurface> ERK;

    
    clock_t begin, end;
    Vector<Real> time1(loop, 0);
    radio = 1;
    center = Point{0.0,0.0};
    n = 128;
    dt = 1e-5;
    Vector<Vector<Point> > list1;
    /*
    cout << "Circle test:" << endl;
    for (int j = 0 ; j < loop ; j++){
      //get the initial curve
      Vector<Point> pts;
      pts.push_back({center[0] + radio, center[1]});
      for (int i = 1; i < n; i++)
        {
          pts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] +
    radio * sin(2 * M_PI / n * i)});
        }
      pts.push_back({center[0] + radio, center[1]});
      cout << "----------- n = " << n << ": ------------" << endl;
      crv = fitCurve<4>(pts, Curve<2, 4>::periodic);
      Vector<Curve<2, 4>> vcrv{crv};
      YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);

      MARS2DIMV<4,VectorOnHypersurface> CM(&ESDIRK4, 4 * M_PI * radio / n,
    test.rtiny);

      VectorForCurvatureFlow<2,4> SFV1;
      VectorForSurfaceDiffusionFlow<2,4> SFV2;


      //Vector<Real> arclength = calArcLength<3>(crv);

      // Vector<Real> LocalFD2coes =  calLocalFDcoes<4>(arclength,3,2);
      // cout << "cal FD2 coes with order 4 , i = 3:" << endl;
      // for (int l = 0 ; l < (int)LocalFD2coes.size() ; l++)
      //   cout << LocalFD2coes[l] << " ";
      // cout << endl;


      // Vector<Point> velocity = SFV2(pts,0);
      // Vector<Point> exactVelocity = CircleCurvatureLaplacian(radio,n);
      // Real maxnorm = MaxNormVelocity(velocity,exactVelocity);
      // cout << "max-norm: " << maxnorm << endl;
      // Real onenorm =
      // OneNormVelocity(velocity,exactVelocity,2*M_PI*radio/n);
      // cout << "1-norm: " << onenorm << endl;


      // list1.push_back(velocity);


      // for (auto it = velocity.begin(); it != velocity.end() ; ++it)
      //   cout << *it << endl;
      // Tensor<Real,2> Jacobi = SFV1.getJacobi(pts,0);
      // for (int l = 0; l < (Jacobi.size())[0]; l++){
      //   for (int m = 0; m < (Jacobi.size())[1]; m++)
      //     cout << Jacobi(l,m) << " ";
      //   cout << endl;
      // }

      begin = clock();
      CM.trackInterface(SFV2, YS, 0, dt, 0.001);
      end = clock();
      time1[j] = (double)(end - begin) / CLOCKS_PER_SEC;
      Curve<2,4> crvn = (YS.getBoundaryCycles())[0];
      crvs.push_back(crvn);
      // Vector<Real> knots = crvn.getKnots();
      // for (int l = 0 ; l < (int)knots.size() ; l++)
      //   cout << crvn(knots[l])[0] << " " << crvn(knots[l])[1] <<
      // endl;

      n*=2;
      dt /= 2;
    }

    //get the approx solution
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

    //output the convergency rate
    auto it1 = crvs.begin();
    //auto it2 = crvs.begin() + loop - 1;
    auto it3 = crvs.end();
    auto result1 = exactError(Vector<Crv>(it1, it3), rcrv, tol);
    //auto result2 = exactError(Vector<Crv>(it1, it2), r2crv, tol);
    //auto result2 = exactError(Vector<Crv>(it2, it3), rcrv, tol);
    cout << "Error and ratio:" << endl;
    for (int l = 0 ; l < (int)result1.size() ; l++)
      cout << result1[l] << endl;
    cout << "CPU time:" << endl;
    for (int l = 0 ; l < (int)time1.size() ; l++)
      cout << time1[l] << endl;

    // n = 64;
    // for (int i = 0; i < loop - 2 ; i++){
    //   cout << "----------- n = " << n << ": ------------" << endl;
    //   cout << "max-norm: " <<
    // MaxNormVelocity(Extrapolation(list1[i],list1[loop-1]),list1[i]) <<
    // endl;
    //   cout << "1-norm: " <<
    OneNormVelocity(Extrapolation(list1[i],list1[loop-1]),list1[i],2*M_PI*radio/n)
    <<
    // endl;
    //   n*=2;
    // }
    */

    // n = 32;
    dt = 1e-4;
    n = 32;
    radio = 0.3;
    center = Point{1.0,1.0};
    Vector<Vector<Point> > list2;
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
      //cout << "----------- n = " << n << ": ------------" << endl;
      crv = fitCurve<4>(pts, Curve<2, 4>::periodic);
      Vector<Curve<2, 4>> vcrv{crv};

      YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);
      MARS2DIMV<4,VectorOnHypersurface> CM(&ESDIRK4, 4* M_PI/ n, test.rtiny);
      
      
      
      
      VectorForCurvatureFlow<2,4> SFV1;
      VectorForSurfaceDiffusionFlow<2,4> SFV2;
      
      // Vector<Point> velocity = SFV1(pts,0);
      // list2.push_back(velocity);
      
      
      // for (auto it = velocity.begin(); it != velocity.end() ; ++it)
      //   cout << *it << endl;
      
      // Tensor<Real,2> Jacobi = SFV2.getJacobi(pts,0);
      // for (int l = 0; l < (Jacobi.size())[0]; l++){
      //   for (int m = 0; m < (Jacobi.size())[1]; m++)
      //     cout << Jacobi(l,m) << " ";
      //   cout << endl;
      // }

      begin = clock();
      CM.trackInterface(SFV2, YS, 0, dt, 6e-2);
      end = clock();
      time1[j] = (double)(end - begin) / CLOCKS_PER_SEC;
      Curve<2,4> crvn = (YS.getBoundaryCycles())[0];
      crvs.push_back(crvn);
      Vector<Real> knots = crvn.getKnots();
      Vector<Point> newpts;
      for (int l = 0 ; l < (int)knots.size() ; l++){
        cout << crvn(knots[l])[0] << " " << crvn(knots[l])[1] <<
      endl;
        newpts.push_back(crvn(knots[l]));
      }
        

      // Tensor<Real,2> Jacobi = SFV2.getJacobi(newpts,0);
      // for (int l = 0; l < (Jacobi.size())[0]; l++){
      //   for (int m = 0; m < (Jacobi.size())[1]; m++)
      //     cout << Jacobi(l,m) << " ";
      //   cout << endl;
      // }
      
      n*=2;
      dt/=2;
    }
    // n = 256;
    // for (int i = 0; i < loop - 2 ; i++){
    //   cout << "----------- n = " << n << ": ------------" << endl;
    //   cout << "max-norm: " <<
    // MaxNormVelocity(Extrapolation(list2[i],list2[loop-1]),list2[i]) <<
    // endl;
    //   cout << "1-norm: " << OneNormVelocity(Extrapolation(list2[i],list2[loop-1]),list2[i],2*M_PI*radio/n) <<
    // endl;
    //   n*=2;
    // }

    //get the approx solution    
    n *= 8; //ensure that the chdlength is smaller than computational solutions'
    Vector<Point> rpts;
    double nradio = sqrt(209.0/200);
    rpts.push_back({nradio, 0.0});
    for (int i = 1; i < n; i++)
    {
        rpts.push_back({nradio * cos(2 * M_PI / n * i), nradio * sin(2 * M_PI / n * i)});
    }
    rpts.push_back({nradio, 0.0});
    auto rcrv = fitCurve<4>(rpts, Curve<2, 4>::periodic);

    //output the convergency rate
    auto it1 = crvs.begin();
    //auto it2 = crvs.begin() + loop - 1;
    auto it3 = crvs.end();
    auto result1 = exactError(Vector<Crv>(it1, it3), rcrv, tol);
    //auto result2 = exactError(Vector<Crv>(it1, it2), r2crv, tol);
    //auto result2 = exactError(Vector<Crv>(it2, it3), rcrv, tol);
    cout << "Error and ratio:" << endl;
    for (int l = 0 ; l < (int)result1.size() ; l++)
      cout << result1[l] << endl;
    for (int l = 0 ; l < (int)time1.size() ; l++)
      cout << time1[l] << endl;
    
    






}

int main()
{
    testFL();

    return 0;
}
