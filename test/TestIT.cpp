#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "InterfaceTracking/MARS2D.h"
#include "InterfaceTracking/ExplicitRK.h"
#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/TestExample.h"

using namespace std;

using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;

int main()
{
    bool plot = false;
    Real tol = 1e-15;
    int stages = 3;
    cout << setiosflags(ios::scientific) << setprecision(4);

    TestIT test = getTest(0);

    //set the initial curve
    int n = test.n;
    Real dt = test.dt;
    int opstride = test.opstride;
    Real radio = test.radio;
    Point center = test.center;
    Vector<Curve<2, 4>> crvs;
    Curve<2, 4> crv;

    ExplicitRungeKutta<2, ClassicRK4> ERK;

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

        if (plot == true)
            CM.trackInterface(*test.velocity, YS, 0, dt, test.T, true, fname, opstride);
        else
            CM.trackInterface(*test.velocity, YS, 0, dt, test.T);

        //get the curve after tracking
        crv = (YS.getBoundaryCycles())[0];
        crvs.push_back(crv);
        n *= 2;
        opstride *= 2;
        dt /= 2;
    }
    
    //get the approx solution
    n *= 8;//ensure that the chdlength is smaller than computational solutions'
    Vector<Point> rpts;
    rpts.push_back({center[0] + radio, center[1]});
    for (int i = 1; i < n; i++)
    {
        rpts.push_back({center[0] + radio * cos(2 * M_PI / n * i), center[1] + radio * sin(2 * M_PI / n * i)});
    }
    rpts.push_back({center[0] + radio, center[1]});
    auto rcrv = fitCurve<4>(rpts, true);

    
    //output the convergency rate
    auto result = exactError(crvs, rcrv, tol);
    //auto result = richardsonError(crvs, tol);
    for (auto &i : result)
    {
        cout << i << "  ";
    }
    cout << endl;

    return 0;
}