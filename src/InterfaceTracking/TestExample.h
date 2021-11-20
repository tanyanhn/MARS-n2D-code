#ifndef _TESTEXAMPLE_
#define _TESTEXAMPLE_

#include "Core/Config.h"
#include "Core/Vec.h"
#include "VectorFunction.h"
#include "VelocityField.h"

using Point = Vec<Real, 2>;

template<class T>
using Vector = std::vector<T>;

class TestIT{
public:
Point cent;
Real radio, dt, rtiny, T;
IT_VectorFunction<2> *velocity;

int n, opstride;
TestIT(Point _cent, Real _radio, int _n, Real _dt, Real _T, IT_VectorFunction<2> *_v, int _opstride, Real _rtiny = 0.1):cent(_cent), radio(_radio), n(_n), dt(_dt), velocity(_v), rtiny(_rtiny), opstride(_opstride), T(_T){};
};

IT_VectorFunction<2> *vortex = new Vortex(8);
IT_VectorFunction<2> *deformation = new Deformation(2);

TestIT getTest(int i)
{
    Vector<TestIT> test;
    //Vortex(8) 
    //n = 64; dt = 0.04; rtiny = 0.01
    //4.26e-5 4.16 2.38e-6 5.04 7.21e-8
    test.push_back(TestIT(Point{0.5, 0.75}, 0.15, 64, 0.04, 8, vortex, 100, 0.01));

    //Deformation(2, 4)
    //n = 128; dt = 0.01; rtiny = 0.01

    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 128, 0.01, 2, deformation, 40, 0.01));

    return test[i];
}  

#endif