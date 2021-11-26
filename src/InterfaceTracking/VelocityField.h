#ifndef _VELOCITYFIELD_H_
#define _VELOCITYFIELD_H_

#include "TimeIntegrator.h"
#include <cmath>

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

class RevRotation : public VectorFunction<2>
{
public:
    RevRotation(Real c1, Real c2, Real _om, Real _T) : rc1(c1), rc2(c2), om(_om), T(_T){};
    const Vec<Real, 2> operator()(const Vec<Real, 2> &pt, Real tn) const
    {
        Vec<Real, 2> y;
        y[0] = om * (rc1 - pt[1]) * cos(M_PI * tn / T);
        y[1] = -om * (rc2 - pt[0]) * cos(M_PI * tn / T);
        return y;
    }

private:
    Real rc1, rc2, om, T;
};

class Vortex : public VectorFunction<2>
{
public:
    Vortex(Real _T) : T(_T){};
    const Vec<Real, 2> operator()(const Vec<Real, 2> &pt, Real tn) const
    {
        Vec<Real, 2> y;
        y[0] = cos(M_PI * tn / T) * pow(sin(M_PI * pt[0]), 2) * sin(2 * M_PI * pt[1]);
        y[1] = -cos(M_PI * tn / T) * pow(sin(M_PI * pt[1]), 2) * sin(2 * M_PI * pt[0]);
        return y;
    }

private:
    Real T;
};

class Deformation : public VectorFunction<2>
{
public:
    Deformation(Real _T, int _n = 4) : T(_T), n(_n){};
    const Vec<Real, 2> operator()(const Vec<Real, 2> &pt, Real tn) const
    {
        Vec<Real, 2> y;
        y[0] = cos(M_PI * tn / T) * sin(n * M_PI * (pt[0] + 0.5)) * sin(n * M_PI * (pt[1] + 0.5));
        y[1] = cos(M_PI * tn / T) * cos(n * M_PI * (pt[0] + 0.5)) * cos(n * M_PI * (pt[1] + 0.5));
        return y;
    }

private:
    Real T;
    int n;
};

#endif