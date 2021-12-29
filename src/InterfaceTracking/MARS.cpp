#include "MARS.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

template <int Dim, int Order>
void MARS<Dim, Order>::trackInterface(const VectorFunction<Dim> &v, YS &ys, Real StartTime, Real dt, Real EndTime, bool output, string fName, int opstride)
{
    Real T = StartTime;
    Real t = dt;
    int step = 1;
    if (output == true)
    {
        ofstream of(string(fName + "_Start.dat"), ios_base::binary);
        ys.dump(of);
    }

    while ((dt > 0 && T < EndTime) || (dt < 0 && T > EndTime))
    {
        if (abs(EndTime - T) < abs(dt))
        {
            t = EndTime - T;
        }
        cout << "Step: " << step << "     timestep: " << t << endl;

        timeStep(v, ys, T, t);

        cout << endl;
        T += t;

        if (output == true && T == EndTime)
        {
            ofstream of(string(fName + "_End.dat"), ios_base::binary);
            ys.dump(of);
            break;
        }
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

template class MARS<2, 4>;
