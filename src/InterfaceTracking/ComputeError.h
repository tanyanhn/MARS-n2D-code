#ifndef _COMPUTEERROR_H_
#define _COMPUTEERROR_H_

#include "XorArea.h"
#include <cmath>
#include <numeric>

using Crv = Curve<2, 4>;
using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;

Vector<Real> circleerror(const Vector<Crv> &crvs, const Point &center, Real radio)
{
    int n = crvs.size();
    Vector<Real> result(2 * n - 1);
    for (int i = 0; i < n; i++)
    {
        auto polys = crvs[i].getPolys();
        auto knots = crvs[i].getKnots();
        int cn = polys.size();
        Vector<Real> dis(cn);
        for (int j = 0; j < cn; j++)
        {
            auto pt = crvs[i]((knots[j] + knots[j + 1]) / 2);
            dis[j] = abs(norm(pt - center, 2) - radio);
        }

        auto sum = accumulate(dis.begin(), dis.end(), 0.0);
        result[2 * i] = sum / (cn + 1);

        /*
        std::cout << std::endl;
        auto maxp = max_element(dis.begin(), dis.end());
        result[2 * i] = *maxp;
        */
    }
    for (int i = 0; i < n - 1; i++)
    {
        result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
    }
    return result;
}

Vector<Real> exactError(const Vector<Crv> &crvs, const Crv &exact, Real tol)
{
    int n = crvs.size();
    Vector<Real> result(2 * n - 1);
    for (int i = 0; i < n; i++)
    {
        result[2 * i] = xorArea(exact, crvs[i], tol);
        //result[2 * i] = xorArea(crvs[i], exact, tol);
    }
    for (int i = 0; i < n - 1; i++)
    {
        result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
    }
    return result;
}

Vector<Real> richardsonError(const Vector<Crv> &crvs, Real tol)
{
    int n = crvs.size();
    Vector<Real> result(2 * n - 3);
    for (int i = 0; i < n - 1; i++)
    {
        //result[2 * i] = xorArea(crvs[i + 1], crvs[i], tol);
        result[2 * i] = xorArea(crvs[i], crvs[i + 1], tol);
    }
    for (int i = 0; i < n - 2; i++)
    {
        result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
    }
    return result;
}

#endif