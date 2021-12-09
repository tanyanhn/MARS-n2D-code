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
    //
    //4.80e-5 4.25 2.52e-6 4.96 8.08e-8 7.41 4.77e-10 4.09 2.80e-11(exact)
    //2.47e-5 3.89 1.66e-6 5.02 5.10e-8 6.80 4.56e-10 4.07 2.72e-11(coarse to fine)
    //time_IML  2.8553e-01  1.9017e+00  1.0669e+00  1.9480e+00  4.1166e+00  1.9841e+00  1.6286e+01  2.0179e+00  6.5956e+01
    //time_vect 2.8354e-01  1.8871e+00  1.0488e+00  1.9537e+00  4.0629e+00  1.9833e+00  1.6064e+01  2.0184e+00  6.5082e+01
    //time_list 3.0822e-01  1.8851e+00  1.1385e+00  1.9478e+00  4.3919e+00  1.9624e+00  1.7116e+01  2.0238e+00  6.9601e+01
    test.push_back(TestIT(Point{0.5, 0.75}, 0.15, 64, 0.04, 8, vortex, "Vortex", 100, 0.01));

    //test 1
    //Deformation(2, 4)
    //n = 128; dt = 0.01; rtiny = 0.01
    //2.44e-5 3.77 1.79e-6 5.23 4.76e-8(richardson)
    //2.74e-5 3.81 1.95e-6 5.23 5.22e-8(exact)
    //time_list 1.73e+1 2.36e+2 3.58e+3
    //time_vect 1.75e+1 2.43e+2 3.77e+3
    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 128, 0.01, 2, deformation, "Deformation", 40, 0.01));

    //test 2
    //Deformation(2, 4)
    //n = 128; dt = 0.01; rtiny = 0.1
    //3.60e-5 4.61 1.47e-6 5.40 3.49e-8(exact)
    //time_list 1.43e+1 1.93e+2 3.99e+3
    //time_vect 1.43e+1 2.01e+2 4.05e+3
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
    //4.93e-5  4.41  2.31e-6  5.10  6.77e-8  6.21  9.12e-10 
    //4.55e-5  4.66  1.80e-6  5.16  5.06e-8  5.85  8.74e-10 (coarse to fine)
    //time_list 3.7478e+00  5.4055e+01  7.8974e+02  1.2461e+04
    //time_vect 3.9465e+00  5.6405e+01  8.6824e+02  1.2688e+04 
    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 128, 0.04, 2, deformation, "Deformation", 40, 0.01));

    //test 6
    //Deformation(2, 4)
    //n = 128; dt = 0.02; rtiny = 0.01
    //2.42e-05  3.66  1.91e-06  5.20  5.19e-08  6.56  5.49e-10  
    //time_list 3.5678e-01  1.3148e+00  5.2866e+00  2.0518e+01
    //time_vect 3.0234e-01  1.2060e+00  4.8464e+00  1.9286e+01
    test.push_back(TestIT(Point{0.5, 0.5}, 0.15, 128, 0.02, 2, deformation, "Deformation", 40, 0.01));

    //test 7
    //Vortex(8)
    //n = 64; dt = 0.02; rtiny = 0.01
    
    test.push_back(TestIT(Point{0.5, 0.75}, 0.15, 64, 0.02, 8, vortex, "Vortex", 100, 0.01));

    return test[i];
}

#endif