#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "InterfaceTracking/MARS.h"
#include "InterfaceTracking/ExplicitRK.h"
#include "InterfaceTracking/TimeIntegrator.h"
#include "InterfaceTracking/VectorFunction.h"

using namespace std;

template <int N>
using rVec = Vec<Real, N>;

template <class T>
using Vector = std::vector<T>;

class V : public IT_VectorFunction<2>
{
public:
    const Vec<Real, 2> operator()(const Vec<Real, 2> &pt, Real TN) const
    {
        Vec<Real, 2> y;
        y[0] = pt[0]/norm(pt);
        y[1] = pt[1]/norm(pt);
        return y;
    }
};

template <template <class T> class cont>
void func(cont<int> i) { return; };

int main()
{
    cout << setiosflags(ios::scientific) << setprecision(4);
    Vector<rVec<2>> pts;
    Real count = 0.0;
    pts.push_back({1, 0});
    //cout << count << endl;
    while (count < 2 * M_PI)
    {
        int num = rand() % 2 + 1;
        count += M_PI / 6 * num / 2;
        if (count < 2 * M_PI)
            pts.push_back({cos(count), sin(count)});
        else
        {
            count = 2 * M_PI;
            pts.push_back({1, 0});
        }
        //cout << count << endl;
    }

    Curve<2,4> crv = fitCurve<4>(pts,true);

    Vector<Curve<2,4>> vcrv{crv};

    YinSet<2,4> YS(SegmentedRealizableSpadjor<4>(vcrv),tol);

    ExplicitRungeKutta<2, RK_ButcherTableau<ERK, ClassicRK4>> ERK;

    rVec<2> pt{1,1};

    MARS<2,4> CM(&ERK, M_PI/6, 0.1);


    cout << "chdLenRange: " << M_PI/6 << ", " << M_PI/6*0.1 << endl;

    auto ptsorigin = pts;
    auto stayid = CM.removeSmallEdges(pts);

    /*CM.trackInterface(V(), YS, 0, 0.1, 1);

    crv = (YS.getBoundaryCycles())[0];

    Vector<Polynomial<4,rVec<2>>> polys = crv.getPolys();

    for(auto i = polys.begin() ; i!=polys.end()-1; i++)
    {
        cout << norm((*(i+1))[0]-(*i)[0]) << endl;
    }
    cout << norm((polys[0])[0]-(*(polys.end()-1))[0]) << endl;*/


    //MARS<2,4> CM(&ERK, )

    return 0;
}