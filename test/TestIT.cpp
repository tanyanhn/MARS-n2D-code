#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "InterfaceTracking/MARS.h"
#include "InterfaceTracking/ExplicitRK.h"
#include "InterfaceTracking/TimeIntegrator.h"
#include "InterfaceTracking/VectorFunction.h"
#include "InterfaceTracking/XorArea.h"
#include "InterfaceTracking/ComputeError.h"
#include "InterfaceTracking/VelocityField.h"

using namespace std;

template <int N>
using rVec = Vec<Real, N>;

template <class T>
using Vector = std::vector<T>;

int main()
{
    cout << setiosflags(ios::scientific) << setprecision(4);

    //set the initial curve
    int n = 16;
    Curve<2, 4> crv;
    Vector<rVec<2>> pts;
    //Real ang = M_PI / n;
    pts.push_back({1, 0});
    for (int i = 1; i < n; i++)
    {
        pts.push_back({cos(2 * M_PI / n * i), sin(2 * M_PI / n * i)});
    }
    pts.push_back({1, 0});
    crv = fitCurve<4>(pts, true);

    //set the CubicMARS method
    Vector<Curve<2, 4>> vcrv{crv};
    YinSet<2, 4> YS(SegmentedRealizableSpadjor<4>(vcrv), tol);
    ExplicitRungeKutta<2, ClassicRK4> ERK;
    MARS<2, 4> CM(&ERK, M_PI / 6, 0.1);

    //start tracking interface
    CM.trackInterface(Deformation(1), YS, 0, 0.01, 1);

    //get the curve after tracking
    crv = (YS.getBoundaryCycles())[0];
    auto polys = crv.getPolys();

    for (auto &i : polys)
    {
        cout << "(" << i(0)[0] << ", " << i(0)[1] << ")" << endl;
    }
    cout << "(" << polys[0](0)[0] << ", " << polys[0](0)[1] << ")" << endl;

    //MARS<2,4> CM(&ERK, )

    return 0;
}