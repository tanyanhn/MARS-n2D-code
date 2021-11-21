#ifndef _COMPUTEERROR_
#define _COMPUTEERROR_

#include "XorArea.h"
#include <cmath>

using Crv = Curve<2, 4>;
using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;

Vector<Real> exactError(const Vector<Crv> &crvs, const Crv &exact, Real tol)
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

/*
Vector<Real> circleError(const Vector<Crv> &crvs, Point center, Real radio)
{
    int n = crvs.size();
    Vector<Real> result(2 * n - 1);

    for (int i = 0; i < n; i++)
    {
        auto polys = crvs[i].getPolys();
        auto knots = crvs[i].getKnots();
        int num = polys.size();

        //compute arclength
        Vector<Real> arcLength(num);
        for (int j = 0; j < num - 1; j++)
        {
            arcLength[j] = radio * acos(dot(polys[j][0] - center, polys[j + 1][0] - center) / norm(polys[j][0] - center, 2) / norm(polys[j + 1][0] - center, 2));
        }
        arcLength[num - 1] = radio * abs(acos(dot(polys[num - 1][0] - center, polys[0][0] - center) / norm(polys[num - 1][0] - center, 2) / norm(polys[0][0] - center)));

        for(auto x: arcLength)
        {
            std::cout << x << " ";
        }
        std::cout << std::endl;

        Vector<Real> normaldist(num);
        for (int j = 0; j < num; j++)
        {
            auto pt = polys[j][0];
            normaldist[j] = norm(pt - center, 2) - radio;
        }

        Real area = 0;
        for (int j = 0; j < num - 1; j++)
        {
            if (normaldist[j] * normaldist[j + 1] >= 0)
                area += std::abs(normaldist[j] + normaldist[j + 1]) / 2 * arcLength[j];
            else
                area += (pow(normaldist[j], 2) + pow(normaldist[j + 1], 2)) * arcLength[j] / 2 / (std::abs(normaldist[j]) + std::abs(normaldist[j + 1]));
        }
        if (normaldist[num - 1] * normaldist[0] >= 0)
            area += std::abs(normaldist[num - 1] + normaldist[0]) / 2 * arcLength[num - 1];
        else
            area += (pow(normaldist[num - 1], 2) + pow(normaldist[0], 2)) * arcLength[num - 1] / 2 / (std::abs(normaldist[num - 1]) + std::abs(normaldist[0]));

        result[2 * i] = area;
    }

    for (int i = 0; i < n - 1; i++)
    {
        result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
    }
    return result;
}
*/

#endif