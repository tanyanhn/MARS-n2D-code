#ifndef _XORAREA_
#define _XORAREA_

#include "Core/Curve.h"
#include "Core/numlib.h"
#include <tuple>
#include <cmath>
#include <numeric>//debug mode 

using Point = Vec<Real, 2>;

using Crv = Curve<2, 4>;

template <int order>
Real localArcLength(const Polynomial<order, Point> &poly, Real start, Real end)
{
    auto dxdt = getComp(poly, 0).der();
    auto dydt = getComp(poly, 1).der();
    auto ndsdt = [&](Real t)
    {
        Vec<Real, 2> ds{dxdt(t), dydt(t)};
        return norm(ds, 2);
    };
    return quad<4>(ndsdt, start, end);
}

std::tuple<bool, Real, Real> segIntSeg(Point lp1, Point hp1, Point lp2, Point hp2, Real tol)
{
    using rVec = Vec<Real, 2>;
    rVec A[2], b;
    A[0] = hp1 - lp1;
    A[1] = hp2 - lp2;
    b = lp2 - lp1;
    Real det = cross(A[0], A[1]);
    //Real sc = norm(A[0]);

    // solve for intersections by Cramer's rule
    Real x[2];
    x[0] = cross(b, A[1]) / det;
    x[1] = cross(A[0], b) / det;

    if (x[0] > -tol && x[0] < 1 + tol && x[1] > -tol && x[1] < 1 + tol)
    {
        return std::make_tuple(true, 2 * (x[0] - 0.5), 2 * (x[1] - 0.5));
    }
    return std::make_tuple(false, 0, 0);
}

Real segIntCurve(Point pt, Point grad, const Crv &crv, Real tol)
{
    auto polys = crv.getPolys();
    auto knots = crv.getKnots();
    int Num = polys.size();
    std::vector<Point> pts(Num + 1);
    for (int i = 0; i < Num; i++)
    {
        pts[i] = polys[i][0];
    }
    pts[Num] = polys[0][0];

    std::vector<std::tuple<int, Real, Real>> intersect;
    for (int i = 0; i < Num; i++)
    {
        std::tuple<bool, Real, Real> ifint = segIntSeg(pt - grad, pt + grad, pts[i], pts[i + 1], tol);
        if (std::get<0>(ifint) == true)
            intersect.push_back(std::make_tuple(i, std::get<1>(ifint), std::get<2>(ifint) * 0.5 + 0.5));
    }
    assert(intersect.size() != 0);

    //find the closest intersect point
    Real min = 1;
    int id;
    int no;
    for (int i = 0; i < (int)intersect.size(); i++)
    {
        if (std::abs(std::get<1>(intersect[i])) <= std::abs(min))
        {
            id = std::get<0>(intersect[i]);
            min = std::get<1>(intersect[i]);
            no = i;
        }
    }

    //find the intersection between segment and local curve, return the signed distance
    Real t0 = (knots[id + 1] - knots[id]) * std::get<2>(intersect[no]);
    if (t0 < (knots[id + 1] - knots[id]) * tol && t0 > (knots[id + 1] - knots[id]) * (1 - tol))
        return min;

    auto xcf = getComp(polys[id], 0);
    auto ycf = getComp(polys[id], 1);
    if (std::abs(grad[0]) < tol * norm(grad, 2))
    {
        xcf[0] -= pt[0];
        Real t = fzero(xcf, xcf.der(), t0, 10, (knots[id + 1] - knots[id]) * tol);
        Point apt = polys[id](t);
        return ((apt[1] - pt[1]) / grad[1] > 0) ? norm(apt[1] - pt[1],2) : -norm(apt[1] - pt[1],2);
    }
    if (std::abs(grad[1]) < tol * norm(grad, 2))
    {
        ycf[0] -= pt[1];
        Real t = fzero(ycf, ycf.der(), t0, 10, (knots[id + 1] - knots[id]) * tol);
        Point apt = polys[id](t);
        return ((apt[0] - pt[0]) / grad[0] > 0) ? norm(apt[0] - pt[0],2) : -norm(apt[0] - pt[0],2);
    }
    Real k = grad[1] / grad[0];
    ycf[0] -= pt[1];
    xcf[0] -= pt[0];
    auto coef = ycf - xcf * k;
    Real t = fzero(coef, coef.der(), t0, 10, (knots[id + 1] - knots[id]) * tol);
    Point apt = polys[id](t);
    if (std::abs(grad[0]) > std::abs(grad[1]))
        return ((apt[0] - pt[0]) / grad[0] > 0) ? norm(apt - pt, 2) : -norm(apt - pt, 2);
    else
        return ((apt[1] - pt[1]) / grad[1] > 0) ? norm(apt - pt, 2) : -norm(apt - pt, 2);
}

Real xorArea(const Curve<2, 4> &fCrv, const Curve<2, 4> &cCrv, Real tol)
{
    auto polys = fCrv.getPolys();
    auto knots = fCrv.getKnots();
    int num = polys.size();

    //compute arclength
    std::vector<Real> arcLength(num);
    for (int i = 0; i < num; i++)
    {
        arcLength[i] = localArcLength(polys[i], 0, knots[i + 1] - knots[i]);
    }

    std::vector<Real> normaldist(num);

    for (int i = 0; i < num; i++)
    {
        auto pt1 = polys[i][0];
        auto dxdt = getComp(polys[i], 0).der();
        auto dydt = getComp(polys[i], 1).der();
        Point grad{dydt[0], -dxdt[0]};
        grad = grad / norm(grad, 2);

        //find approx vertical alignment point
        normaldist[i] = segIntCurve(pt1, grad, cCrv, tol);
    }

    //debug mode
    /*
    for(auto &i: normaldist)
    {
        std::cout << i << ", ";
    }
    std::cout << std::endl;
    */

    Real area = 0;
    for (int i = 0; i < num - 1; i++)
    {
        if (normaldist[i] * normaldist[i + 1] >= 0)
            area += std::abs(normaldist[i] + normaldist[i + 1]) / 2 * arcLength[i];
        else
            area += (pow(normaldist[i], 2) + pow(normaldist[i + 1], 2)) * arcLength[i] / 2 / (std::abs(normaldist[i]) + std::abs(normaldist[i + 1]));
    }
    if (normaldist[num - 1] * normaldist[0] >= 0)
        area += std::abs(normaldist[num - 1] + normaldist[0]) / 2 * arcLength[num - 1];
    else
        area += (pow(normaldist[num - 1], 2) + pow(normaldist[0], 2)) * arcLength[num - 1] / 2 / (std::abs(normaldist[num - 1]) + std::abs(normaldist[0]));
    return area;
}

#endif