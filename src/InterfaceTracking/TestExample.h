#ifndef _TESTEXAMPLE_H_
#define _TESTEXAMPLE_H_

#include "VectorFunction.h"
#include "VelocityField.h"

class TestIT
{
    using Point = Vec<Real, 2>;

    template <class T>
    using Vector = std::vector<T>;

public:
    Point center;//center of initial circle
    Real radio;//radio of initial circle
    Real dt;//time step
    Real rtiny;//rtiny of MARS
    Real T;//total time
    VectorFunction<2> *velocity;//velocity fiel
    int n;//cut the initial circle to n equal parts to generate the initial YinSet
    int opstride;//out put stride, for plot
    std::string name;//file name of the out put target, for plot

    TestIT(Point _center, Real _radio, int _n, Real _dt, Real _T, VectorFunction<2> *_v, std::string _name, int _opstride, Real _rtiny = 0.1) : center(_center), radio(_radio), dt(_dt), rtiny(_rtiny), T(_T), velocity(_v), n(_n), opstride(_opstride), name(_name){};
};

//two velocity fields
VectorFunction<2> *vortex = new Vortex(8);
VectorFunction<2> *deformation = new Deformation(2);

TestIT getTest(int i)
{
    Vector<TestIT> test;

    //test 0
    //Vortex(8)
    //n = 64; dt = 0.04; rtiny = 0.01
    //4.26e-5 4.16 2.38e-6 5.04 7.21e-8(richardson)
    //4.80e-5 4.25 2.52e-6 4.97 8.06e-8(exact)
    test.push_back(TestIT(Point{0.5, 0.75}, 0.15, 64, 0.04, 8, vortex, "Vortex", 100, 0.01));

    //test 1
    //Deformation(2, 4)
    //n = 128; dt = 0.01; rtiny = 0.01
    //2.44e-5 3.77 1.79e-6 5.23 4.76e-8(richardson)
    //2.74e-5 3.81 1.95e-6 5.23 5.22e-8(exact)
    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 128, 0.01, 2, deformation, "Deformation", 40, 0.01));

    //test 2
    //Deformation(2, 4)
    //n = 128; dt = 0.01; rtiny = 0.1
    //3.60e-5 4.61 1.47e-6 5.40 3.49e-8(exact)
    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 128, 0.01, 2, deformation, "Deformation", 40));

    //test 3
    //Deformation(2, 4) plot use
    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 256, 0.005, 2, deformation, "Deformation", 50));

    //test 4
    //Vortex(8) plot use
    test.push_back(TestIT(Point{0.5, 0.75}, 0.15, 128, 0.01, 8, vortex, "Vortex", 100, 0.01));

    //test 5
    //Deformation(2, 4)
    //n = 128; dt = 0.04; rtiny = 0.01
    //5.65e-5 4.61 2.31e-6 5.10 6.75e-8(exact)
    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 128, 0.04, 2, deformation, "Deformation", 40, 0.01));

    return test[i];
}

#endif