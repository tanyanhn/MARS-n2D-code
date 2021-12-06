#ifndef _SOLVEPERIODICTRI_H_
#define _SOLVEPERIODICTRI_H_

#include "Core/Config.h"
#include <vector>
#include <assert.h>

std::vector<Real> solveL(const std::vector<Real> &a, const std::vector<Real> &p, const std::vector<Real> &rhs)
{
    assert((int)a.size() + 1 == (int)p.size() && (int)p.size() == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[0] = rhs[0] / p[0];
    for (int i = 1; i < N; i++)
    {
        result[i] = (rhs[i] - result[i - 1] * a[i - 1]) / p[i];
    }
    return result;
}

std::vector<Real> solveU(const std::vector<Real> &q, const std::vector<Real> &rhs)
{
    assert((int)q.size() + 2 == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[N - 1] = rhs[N - 1];
    result[N - 2] = rhs[N - 2];
    for (int i = N - 3; i >= 0; i--)
    {
        result[i] = rhs[i] - q[i] * result[i + 1];
    }
    return result;
}

std::vector<Real> solveD(const std::vector<Real> &t, Real s, const std::vector<Real> &rhs)
{
    assert((int)t.size() + 1 == (int)rhs.size());
    int N = rhs.size();
    std::vector<Real> result(N);
    result[N - 1] = (rhs[N - 1] - rhs[0] * s) / (1 - s * t[0]);
    for (int i = N - 2; i >= 0; i--)
    {
        result[i] = rhs[i] - t[i] * result[N - 1];
    }
    return result;
}

std::vector<Real> solvePeriodicTri(const std::vector<Real> &a, const std::vector<Real> &b, const std::vector<Real> &c, const std::vector<Real> &rhs)
{
    assert(a.size() == b.size() && b.size() == c.size() && c.size() == rhs.size());
    int N = rhs.size();

    //get a in L
    std::vector<Real> newa;
    newa.assign(a.begin() + 1, a.end());

    //get p, q , r , s in L and Utilde
    std::vector<Real> p(N);
    std::vector<Real> q(N - 2);
    std::vector<Real> r(N - 1);
    p[0] = b[0];
    q[0] = c[0] / p[0];
    r[0] = a[0] / p[0];
    for (int i = 1; i < N - 2; i++)
    {
        p[i] = b[i] - a[i] * q[i - 1];
        q[i] = c[i] / p[i];
        r[i] = -a[i] * r[i - 1] / p[i];
    }
    p[N - 2] = b[N - 2] - a[N - 2] * q[N - 3];
    r[N - 2] = (c[N - 2] - a[N - 2] * r[N - 3]) / p[N - 2];
    p[N - 1] = b[N - 1] - a[N - 1] * r[N - 2];
    Real s = c[N - 1] / p[N - 1];

    //get t in D
    std::vector<Real> t(N - 1);
    t[N - 2] = r[N - 2];
    for (int i = N - 3; i >= 0; i--)
    {
        t[i] = r[i] - q[i] * t[i + 1];
    }

    std::vector<Real> tmp;
    tmp = solveL(newa, p, rhs);
    tmp = solveU(q, tmp);
    return solveD(t, s, tmp);
}

#endif