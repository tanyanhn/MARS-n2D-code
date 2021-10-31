#ifndef _IT_VECTORFUNCTION_
#define _IT_VECTORFUNCTION_

#include "Core/Vec.h"

template <int Dim>
class IT_VectorFunction
{
    using Point = Vec<Real, Dim>;

public:
    virtual ~IT_VectorFunction(){};

    virtual const Point operator()(const Point &pt, Real t) const =0 ;
};

#endif