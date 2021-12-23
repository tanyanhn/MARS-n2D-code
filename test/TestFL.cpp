#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include "GeometricFlow/DiscreteVecCurvatureLaplacian.h"
#include "InterfaceTracking/TestExample.h"

using namespace std;

using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;



void testFL()
{
    cout << setiosflags(ios::scientific) << setprecision(6);

    TestIT test = getTest(6);
    //set the initial curve
    int n;
    int loop = 1;
    //Real dt;
    Real radio = test.radio;
    Point center = test.center;
    Vector<Curve<2, 4>> crvs;
    Curve<2, 4> crv;

    //clock_t begin, end;
 
    n = test.n;
    n /= 4;
    //dt = test.dt;
    
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
      //Vector<Real> arclength2 = calArcLength<3>(crv);
      // Tensor<Real,2> FD2coes2 = calFD2coes<4>(arclength2);
      // cout << "cal FD2 coes with order 4 :" << endl;
      // for (int k = 0 ; k < (int)arclength2.size()-1; k++){
      // for (int l = 0 ; l < 6 ; l++)
      //   cout << FD2coes2(l,k) << " ";
      // cout << endl;
      // }
      
      DiscreteVecCurvatureLaplacian<2,2> SFV1;
      Vector<Point> velocity = SFV1(pts,0);
      // for(auto it = velocity.begin(); it != velocity.end() ; ++it){
      //   cout << *it << endl;
      // }
      
      // cout << "test for Localdads2:" << endl;
      // Tensor<Real,2> Localdads2 = calLocaldads2<4>(pts,crv,1,1e-5);
      // for (int k = 0 ; k < 6; k++){
      //   for (int l = 0 ; l < 5; l++)
      //     cout << Localdads2(k,l) << " ";
      // cout << endl;
      // }
      // cout << "test for Localdads1:" << endl;
      // Tensor<Real,2> Localdads1 = calLocaldads1<6>(pts,crv,8,1e-4);
      // for (int k = 0 ; k < 7; k++){
      //   for (int l = 0 ; l < 6; l++)
      //     cout << Localdads1(k,l) << " ";
      // cout << endl;
      // }
      cout << "test for dads2:" << endl;
      Tensor<Real,3> dads2 = caldads2<4>(pts,crv,1e-5);
      for (int s = 0 ; s < (dads2.size())[0] ; s++){
        for (int k = 0 ; k < 6; k++){
          for (int l = 0 ; l < 5; l++)
            cout << dads2(s,k,l) << " ";
          cout << endl;
        }
        cout << "--------------------------------------" << endl;
      }
      cout << "test for dads1:" << endl;
      Tensor<Real,3> dads1 = caldads1<6>(pts,crv,1e-5);
      for (int s = 0 ; s < (dads1.size())[0] ; s++){
        for (int k = 0 ; k < 7; k++){
          for (int l = 0 ; l < 6; l++)
            cout << dads1(s,k,l) << " ";
          cout << endl;
        }
        cout << "------------------------------------------" << endl;
      }
      
      n*=2;
    }

    // Vector<Real> arc;
    // for (int i = 0 ; i < 8 ; i++)
    //   arc.push_back(i);
    // Tensor<Real,2> FD2coes = calFD2coes<4>(arc);
    // for (int k = 0 ; k < (int)arc.size()-1; k++){
    //   for (int l = 0 ; l < 6 ; l++)
    //     cout << FD2coes(l,k) << " ";
    //   cout << endl;
    // }


   










    
    // Vector<Real> arclength1 = calArcLength<2>(crv);
    // cout << "cal arclength with order 3 :" << endl;
    // for( auto it = arclength1.begin(); it != arclength1.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Real> arclength2 = calArcLength<3>(crv);
    // cout << "cal arclength with order 5 :" << endl;
    // for( auto it = arclength2.begin(); it != arclength2.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Real> arclength3 = calArcLength<4>(crv);
    // cout << "cal arclength with order 7 :" << endl;
    // for( auto it = arclength3.begin(); it != arclength3.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Tensor<Real,2> FD2coes1 = calFD2coes<2>(arclength1);
    // cout << "cal FD2 coes with order 2 :" << endl;
    // for (int i = 0 ; i < (int)arclength1.size()-1; i++){
    //   for (int j = 0 ; j < 4 ; j++)
    //     cout << FD2coes1(j,i) << " ";
    //   cout << endl;
    // }
    // Tensor<Real,2> FD2coes2 = calFD2coes<4>(arclength2);
    // cout << "cal FD2 coes with order 4 :" << endl;
    // for (int i = 0 ; i < (int)arclength2.size()-1; i++){
    //   for (int j = 0 ; j < 6 ; j++)
    //     cout << FD2coes2(j,i) << " ";
    //   cout << endl;
    // }
    // Tensor<Real,2> FD2coes3 = calFD2coes<6>(arclength3);
    // cout << "cal FD2 coes with order 6 :" << endl;
    // for (int i = 0 ; i < (int)arclength3.size()-1; i++){
    //   for (int j = 0 ; j < 8 ; j++)
    //     cout << FD2coes3(j,i) << " ";
    //   cout << endl;
    // }
    // Tensor<Real,2> FD1coes1 = calFD1coes<2>(arclength1);
    // cout << "cal FD1 coes with order 2 :" << endl;
    // for (int i = 0 ; i < (int)arclength1.size()-1; i++){
    //   for (int j = 0 ; j < 3 ; j++)
    //     cout << FD1coes1(j,i) << " ";
    //   cout << endl;
    // }
    // Tensor<Real,2> FD1coes2 = calFD1coes<4>(arclength2);
    // cout << "cal FD1 coes with order 4 :" << endl;
    // for (int i = 0 ; i < (int)arclength2.size()-1; i++){
    //   for (int j = 0 ; j < 5 ; j++)
    //     cout << FD1coes2(j,i) << " ";
    //   cout << endl;
    // }
    // Tensor<Real,2> FD1coes3 = calFD1coes<6>(arclength3);
    // cout << "cal FD1 coes with order 6 :" << endl;
    // for (int i = 0 ; i < (int)arclength3.size()-1; i++){
    //   for (int j = 0 ; j < 7 ; j++)
    //     cout << FD1coes3(j,i) << " ";
    //   cout << endl;
    // }
    // Vector<Point> der11 = calDer<2>(pts,crv,1);
    // cout << "cal der1 with order 2 :" << endl;
    // for( auto it = der11.begin(); it != der11.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Point> der12 = calDer<4>(pts,crv,1);
    // cout << "cal der1 with order 4 :" << endl;
    // for( auto it = der12.begin(); it != der12.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Point> der13 = calDer<6>(pts,crv,1);
    // cout << "cal der1 with order 6 :" << endl;
    // for( auto it = der13.begin(); it != der13.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Point> der21 = calDer<2>(pts,crv,2);
    // cout << "cal der2 with order 2 :" << endl;
    // for( auto it = der21.begin(); it != der21.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Point> der22 = calDer<4>(pts,crv,2);
    // cout << "cal der2 with order 4 :" << endl;
    // for( auto it = der22.begin(); it != der22.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Point> der23 = calDer<6>(pts,crv,2);
    // cout << "cal der2 with order 6 :" << endl;
    // for( auto it = der23.begin(); it != der23.end() ; ++it){
    //   cout << *it << endl;
    // }
    // const int num2 = der22.size();
    // Vector<Real> kappa1(num2);
    // for (int i = 0 ; i < num2 ; i++)
    //   kappa1[i] = -der22[i][0]*der12[i][1]+der22[i][1]*der12[i][0];
    // Vector<Real> d2kappa1 = calDer<4>(kappa1,crv,2);
    // cout << "cal d2kappa with order 2:" << endl;
    // for( auto it = d2kappa1.begin(); it != d2kappa1.end() ; ++it){
    //   cout << *it << endl;
    // }
    // Vector<Real> kappa2(num2);
    // for (int i = 0 ; i < num2 ; i++)
    //   kappa2[i] = -der23[i][0]*der13[i][1]+der23[i][1]*der13[i][0];
    // Vector<Real> d2kappa2 = calDer<6>(kappa2,crv,2);
    // cout << "cal d2kappa with order 4:" << endl;
    // for( auto it = d2kappa2.begin(); it != d2kappa2.end() ; ++it){
    //   cout << *it << endl;
    // }
}

int main()
{
    testFL();

    return 0;
}
