#ifndef _COMPUTEERROR_
#define _COMPUTEERROR_

#include "XorArea.h"

using Crv = Curve<2, 4>;

template <class T>
using Vector = std::vector<T>;

Vector<Real> exatError(const Vector<Crv> &crvs, const Crv &exact, Real tol)
{
    int n = crvs.size();
    Vector<Real> result(2 * n - 1);
    for (int i = 0; i < n; i++)
    {
        result[2 * i] = xorArea(exact, crvs[i], tol);
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
        result[2 * i] = xorArea(crvs[i + 1], crvs[i], tol);
    }
    for (int i = 0; i < n - 2; i++)
    {
        result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
    }
    return result;
}

#endif