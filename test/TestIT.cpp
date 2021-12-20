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
#include "InterfaceTracking/MARS2DIMV.h"
#include "InterfaceTracking/ERK.h"
#include "InterfaceTracking/DIRK.h"
#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/TestExample.h"

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
    int loop = 1;
    bool plot = false;
    Real tol = 1e-15;
    int stages = 5;
    cout << setiosflags(ios::scientific) << setprecision(4);

    TestIT test = getTest(6);

    //set the initial curve
    int n;
    Real dt;
    int opstride;
    Real radio = test.radio;
    Point center = test.center;
    Vector<Curve<2, 4>> crvs;
    Curve<2, 4> crv;

    ERK<2, RK::ClassicRK4> ERK;

    Vector<Real> time1(2 * stages - 1, 0);
    Vector<Real> time2(2 * stages - 1, 0);
    clock_t begin, end;
    for (int lp = 0; lp < loop; lp++)
    {
        
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
            MARS2DIMV<4> CM(&ERK, 4 * M_PI * radio / n, test.rtiny);

            ostringstream tmps;
            tmps << k;
            string fname = "results" + test.name + "/No" + tmps.str();

            begin = clock();
            if (plot == true)
                CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname, opstride);
            else
                CM.trackInterface(*test.velocity, YS, 0, dt, test.T);
            end = clock();
            time1[2 * k] += (double)(end - begin) / CLOCKS_PER_SEC;
            //get the curve after tracking
            crv = (YS.getBoundaryCycles())[0];
            crvs.push_back(crv);
            n *= 2;
            opstride *= 2;
            dt /= 2;
        }
        
        //crvs[0] = output(crv, center, radio);
        /*
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
            MARS2D<4> CM(&ERK, 4 * M_PI * radio / n, test.rtiny);

            ostringstream tmps;
            tmps << k;
            string fname = "results" + test.name + "/No" + tmps.str();

            begin = clock();
            if (plot == true)
                CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname, opstride);
            else
                CM.trackInterface(*test.velocity, YS, 0, dt, test.T);
            end = clock();
            time2[2 * k] += (double)(end - begin) / CLOCKS_PER_SEC;
            //get the curve after tracking
            crv = (YS.getBoundaryCycles())[0];
            crvs.push_back(crv);
            n *= 2;
            opstride *= 2;
            dt /= 2;
        }
        */
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

    //output the convergency rate
    auto it1 = crvs.begin();
    auto it2 = crvs.begin() + stages;
    //auto it3 = crvs.end();
    auto result1 = exactError(Vector<Crv>(it1, it2), rcrv, tol);
    //auto result2 = exactError(Vector<Crv>(it1, it2), r2crv, tol);
    //auto result2 = exactError(Vector<Crv>(it2, it3), rcrv, tol);

    //auto result1 = richardsonError(crvs, tol);
    for (auto &i : result1)
    {
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

    for (auto &i : time1)
    {
        i = i / loop;
    }

    for (auto &i : time2)
    {
        i = i / loop;
    }

    for (int i = 0; i < stages - 1; i++)
    {
        time1[2 * i + 1] = log(time1[2 * i + 2]/time1[2 * i]) / log(2);
        time2[2 * i + 1] = log(time2[2 * i + 2]/time2[2 * i]) / log(2);
    }
    cout << "method1 time: ";
    for (auto &i : time1)
    {
        cout << i << "  ";
    }
    cout << endl;

    cout << "method2 time: ";
    for (auto &i : time2)
    {
        cout << i << "  ";
    }
    cout << endl;
}

int main()
{
    testIT();

    return 0;
}