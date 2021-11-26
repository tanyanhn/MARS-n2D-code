#ifndef _VECTORFUNCTION_H_
#define _VECTORFUNCTION_H_

#include "Core/Config.h"
#include "Core/Vec.h"

template <int Dim>
class VectorFunction
{
    using Point = Vec<Real, Dim>;

public:
    virtual ~VectorFunction(){};

    virtual const Point operator()(const Point &pt, Real t) const =0 ;
};

#endif