#ifndef _VECTORFUNCTION_H_
#define _VECTORFUNCTION_H_

#include "Core/Config.h"
#include "Core/Vec.h"
#include "Core/Tensor.h"
#include <vector>

template <int Dim>
class VectorFunction
{
private:
    using Point = Vec<Real, Dim>;

    template <class T>
    using Vector = std::vector<T>;

    virtual const Point operator()(Point pt, Real t) const = 0;

    virtual const Vector<Real> getJacobi(Point pt, Real t) const
    {
        return Vector<Real>(Dim * Dim);
    }

public:
    virtual ~VectorFunction(){};

    virtual const Vector<Point> operator()(const Vector<Point> &pts, Real t = 0) const
    {
        int num = pts.size();
        Vector<Point> vel(num);
        for (int i = 0; i < num; i++)
        {
            vel[i] = this->operator()(pts[i], t);
        }
        return vel;
    }

    virtual const Tensor<Real, 2> getJacobi(const Vector<Point> &pts, Real t = 0) const
    {
        int num = pts.size() - 1;
        Tensor<Real, 2> jacobi(Dim * num);
        Vector<Real> ptjac;

        for (int i = 0; i < num; i++)
        {
            ptjac = getJacobi(pts[i], t);
            for (int j = 0; j < Dim; j++)
            {
                for (int k = 0; k < Dim; k++)
                {
                    jacobi(i + j * num, i + k * num) = ptjac[j * Dim + k];
                }
            }
        }
        return jacobi;
    }
};

#endif