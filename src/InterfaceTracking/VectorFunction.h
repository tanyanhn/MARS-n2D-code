#ifndef _IT_VECTORFUNCTION_
#define _IT_VECTORFUNCTION_

#include "Core/Vec.h"

template <int N>
class IT_VectorFunction
{
    using rVec = Vec<Real, N>;

public:
    virtual ~IT_VectorFunction(){};

    virtual const rVec operator()(const rVec &pt, Real t) const =0 ;
};

#endif