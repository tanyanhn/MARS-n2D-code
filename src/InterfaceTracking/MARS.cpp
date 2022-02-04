#include "MARS.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

template <int Dim, int Order, template <int> class VelocityField>
void MARS<Dim, Order, VelocityField>::trackInterface(const VelocityField<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime, bool output, string fName, int opstride)
{
    Real T = StartTime;
    int stages = ceil(abs(EndTime - StartTime) / abs(dt));
    Real k = (EndTime - StartTime) / stages;
    int step = 1;
    if (output == true)
    {
        ofstream of(string(fName + "_Start.dat"), ios_base::binary);
        ys.dump(of);
    }

    while (step <= stages)
    {
        cout << "Step: " << step << "     timestep: " << k << endl;

        timeStep(v, ys, T, k);

        //
        auto crv = (ys.getBoundaryCycles())[0];
        auto knots = crv.getKnots();
        auto polys = crv.getPolys();
        int ln = polys.size();

        for (int i = 0; i < ln; i++)
        {
            cout << crv(knots[i]) << "  " << crv((knots[i] + knots[i + 1]) / 2) << "  ";
        }
        cout << endl;
        //

        cout << endl;
        T += k;

        if (output == true && step % opstride == 0)
        {
            ostringstream tmps;
            tmps << step;
            ofstream of(string(fName + "_Step" + tmps.str() + ".dat"), ios_base::binary);
            ys.dump(of);
        }
        step++;
    }
    return;
}

template class MARS<2, 2, VectorFunction>;
template class MARS<2, 4, VectorFunction>;

template class MARS<2, 2, VectorOnHypersurface>;
template class MARS<2, 4, VectorOnHypersurface>;
