#ifndef _VELOCITYFIELD_H_
#define _VELOCITYFIELD_H_

#include "TimeIntegrator.h"
#include <cmath>
/*
class Translation : public VectorFunction<2>
{
public:
    Translation(Real u1, Real u2) : v1(u1), v2(u2){};
    const Vec<Real, 2> operator()(const Vec<Real, 2> &pt, Real tn) const
    {
        Vec<Real, 2> y;
        y[0] = v1;
        y[1] = v2;
        return y;
    }

private:
    Real v1, v2;
};
*/
/*
class Rotation : public VectorFunction<2>
{
public:
    Rotation(Real c1, Real c2, Real _om) : rc1(c1), rc2(c2), om(_om){};
    const Vec<Real, 2> operator()(const Vec<Real, 2> &pt, Real tn) const
    {
        Vec<Real, 2> y;
        y[0] = om * (rc1 - pt[1]);
        y[1] = -om * (rc2 - pt[0]);
        return y;
    }

private:
    Real rc1, rc2, om;
};
*/

class RevRotation : public VectorFunction<2>
{
private:
    using Point = Vec<Real, 2>;

    template <class T>
    using Vector = std::vector<T>;

    const Point operator()(const Point &pt, Real t) const
    {
        Point y;
        y[0] = om * (rc1 - pt[1]) * cos(M_PI * t / T);
        y[1] = -om * (rc2 - pt[0]) * cos(M_PI * t / T);
        return y;
    }

    const Vector<Real> getJacobi(Point pt, Real t) const
    {
        Vector<Real> ptjac(4);
        ptjac[0] = 0;
        ptjac[1] = -om * cos(M_PI * t / T);
        ptjac[2] = om * cos(M_PI * t / T);
        ptjac[3] = 0;
        return ptjac;
    }

public:
    RevRotation(Real c1, Real c2, Real _om, Real _T) : rc1(c1), rc2(c2), om(_om), T(_T){};

private:
    Real rc1, rc2, om, T;
};

class Vortex : public VectorFunction<2>
{
private:
    using Point = Vec<Real, 2>;

    template <class T>
    using Vector = std::vector<T>;

    const Point operator()(Point pt, Real t) const
    {
        Point y;
        y[0] = cos(M_PI * t / T) * pow(sin(M_PI * pt[0]), 2) * sin(2 * M_PI * pt[1]);
        y[1] = -cos(M_PI * t / T) * pow(sin(M_PI * pt[1]), 2) * sin(2 * M_PI * pt[0]);
        return y;
    }

    const Vector<Real> getJacobi(Point pt, Real t) const
    {
        Vector<Real> ptjac(4);
        ptjac[0] = cos(M_PI * t / T) * 2 * M_PI * sin(M_PI * pt[0]) * cos(M_PI * pt[0]) * sin(2 * M_PI * pt[1]);
        ptjac[1] = cos(M_PI * t / T) * 2 * M_PI * pow(sin(M_PI * pt[0]), 2) * cos(2 * M_PI * pt[1]);
        ptjac[2] = -cos(M_PI * t / T) * 2 * M_PI * sin(M_PI * pt[1]) * cos(M_PI * pt[1]) * sin(2 * M_PI * pt[0]);
        ptjac[3] = -cos(M_PI * t / T) * 2 * M_PI * pow(sin(M_PI * pt[1]), 2) * cos(2 * M_PI * pt[0]);
        return ptjac;
    }

public:
    Vortex(Real _T) : T(_T){};

private:
    Real T;
};

class Deformation : public VectorFunction<2>
{
private:
    using Point = Vec<Real, 2>;

    template <class T>
    using Vector = std::vector<T>;

    const Point operator()(Point pt, Real t) const
    {
        Point y;
        y[0] = cos(M_PI * t / T) * sin(n * M_PI * (pt[0] + 0.5)) * sin(n * M_PI * (pt[1] + 0.5));
        y[1] = cos(M_PI * t / T) * cos(n * M_PI * (pt[0] + 0.5)) * cos(n * M_PI * (pt[1] + 0.5));
        return y;
    }

    const Vector<Real> getJacobi(Point pt, Real t) const
    {
        Vector<Real> ptjac(4);
        ptjac[0] = cos(M_PI * t / T) * n * M_PI * cos(n * M_PI * (pt[0] + 0.5)) * sin(n * M_PI * (pt[1] + 0.5));
        ptjac[1] = cos(M_PI * t / T) * n * M_PI * sin(n * M_PI * (pt[0] + 0.5)) * cos(n * M_PI * (pt[1] + 0.5));
        ptjac[2] = -cos(M_PI * t / T) * n * M_PI * cos(n * M_PI * (pt[0] + 0.5)) * sin(n * M_PI * (pt[1] + 0.5));
        ptjac[3] = -cos(M_PI * t / T) * n * M_PI * sin(n * M_PI * (pt[0] + 0.5)) * cos(n * M_PI * (pt[1] + 0.5));
        return ptjac;
    }

public:
    Deformation(Real _T, int _n = 4) : T(_T), n(_n){};

private:
    Real T;
    int n;
};

#endif