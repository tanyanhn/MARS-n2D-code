#ifndef _IT_TIMEINTEGRATOR_
#define _IT_TIMEINTEGRATOR_
#include <vector>
#include "Core/Config.h"
#include "Core/Vec.h"
#include "VectorFunction.h"

template<class T>
struct identity
{
    using type = T;
};

template <int Dim>
class IT_TimeIntegrator
{
    template <class T>
    using Vector = std::vector<T>;

    using Point = Vec<Real, Dim>;

public:
    virtual ~IT_TimeIntegrator() {}

    virtual Point timeStep(const IT_VectorFunction<Dim> &v, const Point &pt, Real tn, Real dt) = 0;

    template <class Container>
    void multyTimeStep(const IT_VectorFunction<Dim> &v, Container &pts, Real tn, Real dt)
    {
        typeTimeStep(identity<Container>(), v, pts, tn, dt);
    }

    template<class T>
    void typeTimeStep(identity<T>, const IT_VectorFunction<Dim> &v, T &pts, Real tn, Real dt);

    void typeTimeStep(identity<Vector<Point>>, const IT_VectorFunction<Dim> &v, Vector<Point> &pts, Real tn, Real dt)
    {
        for (auto &i: pts)
        {
            i = timeStep(v, i, tn, dt);
        }
        return;
    }
};

#endif