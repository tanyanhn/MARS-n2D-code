#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "InterfaceTracking/MARS.h"
#include "InterfaceTracking/ExplicitRK.h"
#include "InterfaceTracking/TimeIntegrator.h"
#include "InterfaceTracking/VectorFunction.h"
#include "InterfaceTracking/XorArea.h"
#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/VelocityField.h"
#include "YinSet/YinSet.h"
#include "InterfaceTracking/TestExample.h"

using namespace std;

template <int N>
using rVec = Vec<Real, N>;

template <class T>
using Vector = std::vector<T>;

int main()
{
    Real tol = 1e-15;
    cout << setiosflags(ios::scientific) << setprecision(4);

    TestIT test = getTest(2);

    //set the initial curve
    int n = test.n;
    Real dt = test.dt;
    int opstride = test.dt;
    Real radio = test.radio;
    Point center = test.cent;
    Vector<Curve<2, 4>> crvs;
    Curve<2, 4> crv;

    for (int k = 0; k < 3; k++)
    {
        Vector<rVec<2>> pts;
        //Real ang = M_PI / n;
        pts.push_back({center[0] + radio, center[1]});
        for (int i = 1; i < n; i++)
        {
            pts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
        }
        pts.push_back({center[0] + radio, center[1]});
        crv = fitCurve<4>(pts, true);

        //set the CubicMARS method
        Vector<Curve<2, 4>> vcrv{crv};
        YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);
        ExplicitRungeKutta<2, ClassicRK4> ERK;
        MARS<2, 4> CM(&ERK, 4 * M_PI * radio / n, test.rtiny);

        ostringstream tmps;
        tmps << k;
        //string fname = "resultsVortex/No" + tmps.str();
        string fname = "results" + test.name + "/No" + tmps.str();

        //start tracking interface
        //CM.trackInterface(Rotation(-2, 0, 2 * M_PI), YS, 0, dt, 1);
        //CM.trackInterface(RevRotation(-2, 0, 2 * M_PI, 1), YS, 0, dt, 1);
        //CM.trackInterface(Vortex(8), YS, 0, dt, 8, true, fname, opstride);
        //CM.trackInterface(Deformation(2), YS, 0, dt, 2, true, fname, opstride);

        //CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname, opstride);
        CM.trackInterface(*test.velocity, YS, 0, dt, test.T);

        //get the curve after tracking
        crv = (YS.getBoundaryCycles())[0];
        crvs.push_back(crv);
        n *= 2;
        opstride *= 2;
        dt /= 2;
    }
    //output the convergency rate

    n *= 8;
    Vector<rVec<2>> rpts;
    //Real ang = M_PI / n;
    rpts.push_back({center[0] + radio, center[1]});
    for (int i = 1; i < n; i++)
    {
        rpts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
    }
    rpts.push_back({center[0] + radio, center[1]});
    auto rcrv = fitCurve<4>(rpts, true);

    //auto result = richardsonError(crvs, tol);
    auto result = exactError(crvs, rcrv, tol);
    for (auto &i : result)
    {
        cout << i << "  ";
    }
    cout << endl;

    //output the points
    /*
    for (auto &crv : crvs)
    {
        auto polys = crv.getPolys();

        for (auto &i : polys)
        {
            cout << "(" << i(0)[0] << ", " << i(0)[1] << ") ";
        }
        cout << "(" << polys[0](0)[0] << ", " << polys[0](0)[1] << ")" << endl
             << endl;
    }
    */
    return 0;
}