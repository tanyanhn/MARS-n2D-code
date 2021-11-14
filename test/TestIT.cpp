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
        y[0] = pt[0] / norm(pt);
        y[1] = pt[1] / norm(pt);
        return y;
    }
};

template <template <class T> class cont>
void func(cont<int> i) { return; };

int main()
{
    cout << setiosflags(ios::scientific) << setprecision(4);
    int n = 16;
    Vector<Curve<2,4>> crvs(4);
    for (int k = 0; k < 4; k++)
    {
        Vector<rVec<2>> pts;
        Vector<rVec<2>> pts2;
        Real ang = M_PI / n;
        pts.push_back({1, 0});
        pts2.push_back({cos(ang), sin(ang)});
        for (int i = 1; i < n; i++)
        {
            pts.push_back({cos(2 * M_PI / n * i), sin(2 * M_PI / n * i)});
            pts2.push_back({cos(2 * M_PI / n * i + ang), sin(2 * M_PI / n * i + ang)});
        }
        pts.push_back({1, 0});
        pts2.push_back({cos(ang), sin(ang)});

        crvs[k] = fitCurve<4>(pts, true);

        n *= 2;
    }

    auto result = richardsonError(crvs, 1e-10);

    for(auto &i: result)
    {
        cout << i << " ";
    }
    cout << endl;

    /*Vector<Curve<2,4>> vcrv{crv};

    YinSet<2,4> YS(SegmentedRealizableSpadjor<4>(vcrv),tol);

    ExplicitRungeKutta<2, ClassicRK4> ERK;

    rVec<2> pt{1,1};

    MARS<2,4> CM(&ERK, M_PI/6, 0.1);


    cout << "chdLenRange: " << M_PI/6 << ", " << M_PI/6*0.1 << endl;


    CM.trackInterface(V(), YS, 0, 0.1, 1);

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