#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include "InterfaceTracking/MARS2D.h"
#include "InterfaceTracking/MARS2DLIST.h"
#include "InterfaceTracking/ExplicitRK.h"
#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/TestExample.h"
#include "InterfaceTracking/SolvePeriodicTri.h"

using namespace std;

using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;

Crv output(const Crv &crv, Point center, Real radio)
{
    Vector<Point> npts;
    auto polys = crv.getPolys();
    auto knots = crv.getKnots();
    int m = 32;
    for (int i = 0; i < (int)polys.size(); i++)
    {
        for (int j = 0; j < m; j++)
        {
            Real dt = j * (knots[i + 1] - knots[i]) / m;
            auto pt = polys[i](dt);
            npts.push_back(pt);
        }
    }
    npts.push_back(polys[0][0]);
    return fitCurve<4>(npts, true);
}

void testIT()
{
    bool plot = false;
    Real tol = 1e-15;
    int stages = 4;
    cout << setiosflags(ios::scientific) << setprecision(4);

    TestIT test = getTest(6);

    //set the initial curve
    int n = test.n;
    Real dt = test.dt;
    int opstride = test.opstride;
    Real radio = test.radio;
    Point center = test.center;
    Vector<Curve<2, 4>> crvs;
    Curve<2, 4> crv;

    ExplicitRungeKutta<2, ClassicRK4> ERK;

    Vector<Real> time1(stages);
    Vector<Real> time2(stages);
    clock_t begin, end;

    for (int k = 0; k < stages; k++)
    {
        //get the initial curve
        Vector<Point> pts;
        pts.push_back({center[0] + radio, center[1]});
        for (int i = 1; i < n; i++)
        {
            pts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
        }
        pts.push_back({center[0] + radio, center[1]});
        crv = fitCurve<4>(pts, true);
        Vector<Curve<2, 4>> vcrv{crv};
        YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);

        //set the CubicMARS method
        MARS2D<4, list> CM(&ERK, 4 * M_PI * radio / n, test.rtiny);

        ostringstream tmps;
        tmps << k;
        string fname = "results" + test.name + "/No" + tmps.str();

        begin = clock();
        if (plot == true)
            CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname, opstride);
        else
            CM.trackInterface(*test.velocity, YS, 0, dt, test.T);
        end = clock();
        time1[k] = (double)(end - begin) / CLOCKS_PER_SEC;
        //get the curve after tracking
        crv = (YS.getBoundaryCycles())[0];
        crvs.push_back(crv);
        n *= 2;
        opstride *= 2;
        dt /= 2;
    }

    //crvs[0] = output(crv, center, radio);

    
    n = test.n;
    dt = test.dt;
    opstride = test.opstride;
    

    for (int k = 0; k < stages; k++)
    {
        //get the initial curve
        Vector<Point> pts;
        pts.push_back({center[0] + radio, center[1]});
        for (int i = 1; i < n; i++)
        {
            pts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
        }
        pts.push_back({center[0] + radio, center[1]});
        crv = fitCurve<4>(pts, true);
        Vector<Curve<2, 4>> vcrv{crv};
        YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);

        //set the CubicMARS method
        MARS2D<4, vector> CM(&ERK, 4 * M_PI * radio / n, test.rtiny);

        ostringstream tmps;
        tmps << k;
        string fname = "results" + test.name + "/No" + tmps.str();

        begin = clock();
        if (plot == true)
            CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname, opstride);
        else
            CM.trackInterface(*test.velocity, YS, 0, dt, test.T);
        end = clock();
        time2[k] = (double)(end - begin) / CLOCKS_PER_SEC;
        //get the curve after tracking
        crv = (YS.getBoundaryCycles())[0];
        crvs.push_back(crv);
        n *= 2;
        opstride *= 2;
        dt /= 2;
    }
    

    //get the approx solution

    n *= 8; //ensure that the chdlength is smaller than computational solutions'
    Vector<Point> rpts;
    rpts.push_back({center[0] + radio, center[1]});
    for (int i = 1; i < n; i++)
    {
        rpts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
    }
    rpts.push_back({center[0] + radio, center[1]});
    auto rcrv = fitCurve<4>(rpts, true);
    /*
    n *= 2; //ensure that the chdlength is smaller than computational solutions'
    Vector<Point> r2pts;
    r2pts.push_back({center[0] + radio, center[1]});
    for (int i = 1; i < n; i++)
    {
        r2pts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
    }
    r2pts.push_back({center[0] + radio, center[1]});
    auto r2crv = fitCurve<4>(r2pts, true);
    */
    //output the convergency rate
    auto it1 = crvs.begin();
    auto it2 = crvs.begin() + stages;
    auto it3 = crvs.end();
    auto result1 = exactError(Vector<Crv>(it1, it2), rcrv, tol);
    //auto result2 = exactError(Vector<Crv>(it1, it2), r2crv, tol);
    auto result2 = exactError(Vector<Crv>(it2, it3), rcrv, tol);

    //auto result1 = richardsonError(crvs, tol);
    for (auto &i : result1)
    {
        cout << i << "  ";
    }
    cout << endl;
    
    for (auto &i : result2)
    {
        cout << i << "  ";
    }
    cout << endl;
    
    
    cout << "  list time: ";
    for (auto &i:time1)
    {
        cout << i << "  ";
    }
    cout << endl;

    cout << "vector time: ";
    for (auto &i:time2)
    {
        cout << i << "  ";
    }
    cout << endl;
}

int main()
{
    //testIT();
    Vector<Real> a{1,2,3,4};
    Vector<Real> b{11,12,13,14};
    Vector<Real> c{21,22,23,24};
    Vector<Real> rhs{1,4,7,10};

    Vector<Real> res = solvePeriodicTri(a,b,c,rhs);
    cout << b[0]*res[0] + c[0]*res[1] + a[0]*res[3] << endl;
    cout << a[1]*res[0] + b[1]*res[1] + c[1]*res[2] << endl;
    cout << a[2]*res[1] + b[2]*res[2] + c[2]*res[3] << endl;
    cout << a[3]*res[2] + b[3]*res[3] + c[3]*res[0] << endl;


    return 0;
}