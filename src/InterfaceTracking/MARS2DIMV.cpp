#include "MARS2DIMV.h"
#include "Core/Polynomial.h"
#include <fstream>
#include <sstream>
#include <list>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <type_traits>

using namespace std;

using Point = Vec<Real, 2>;

template <class T>
using Vector = vector<T>;

static Real tol = 1e-15;

template <int Order>
void MARS2DIMV<Order>::discreteFlowMap(const VectorFunction<2> &v, Vector<Point> &pts, Real tn, Real dt)
{
    Base::TI->timeStep(v, pts, tn, dt);
    return;
}

/**
 * @brief helper function of removeSmallEdges
 * @tparam Container 
 * @param  ids              id of the current points
 * @param  pts              
 * @param  lowBound         
 * @return has removed points or not
 */
bool removeIMV(Vector<unsigned int> &ids, Vector<Point> &pts, Real lowBound)
{
    //set distances
    int num = pts.size();

    int count = 1;
    int i = 1;

    //function<void(int&, int&)> addpi;

    //lambda: split between two points
    auto addpi = [&]() -> void
    {
        if (count != i)
        {
            pts[count] = pts[i];
            ids[count] = ids[i];
        }
        count++;
        i++;
        return;
    };

    Vector<Real> dist(num - 1);
    for (int j = 0; j < num - 1; j++)
    {
        dist[j] = norm(pts[j + 1] - pts[j], 2);
    }

    bool predelete = false; //mark whether the pre-point has been erased
    if (dist[0] < lowBound)
    {
        predelete = true;
        i++;
    }

    while (i < num - 2)
    {
        if (dist[i] >= lowBound)
        {
            predelete = false;
            addpi();
        }
        else
        {
            if (dist[i - 1] <= dist[i + 1])
            {
                if (predelete == false)
                {
                    i++;
                    predelete = true;
                }
                else
                {
                    addpi();
                    predelete = false;
                }
            }
            else
            {
                addpi();
                i++;
                predelete = true;
            }
        }
    }
    if (i == num - 1)
    {
        addpi();
        pts.resize(count);
        ids.resize(count);
        return num != count;
    }
    else
    {
        if (predelete == true || dist[i] >= lowBound)
        {
            addpi();
            addpi();
            pts.resize(count);
            ids.resize(count);
            return num != count;
        }
        i++;
        addpi();
        pts.resize(count);
        ids.resize(count);
        return true;
    }
}

template <int Order>
Vector<unsigned int> MARS2DIMV<Order>::removeSmallEdges(Vector<Point> &pts)
{
    int num = pts.size();
    Vector<unsigned int> ids(num);
    for (int i = 0; i < num; i++)
    {
        ids[i] = i;
    }

    while (true)
    {
        //output ids, debug use
        std::cout << "(";
        for (auto &id : ids)
        {
            std::cout << id << ", ";
        }
        std::cout << "\b)" << std::endl;

        if (!removeIMV(ids, pts, (chdLenRange.lo())[0]))
            break;
    }

    Vector<unsigned int> rids(num - ids.size());
    int it = 0;
    unsigned int allid = 0;
    int count = 0;
    while (count < (int)rids.size() - 1)
    {
        if (ids[it] != allid)
        {
            rids[count] = allid;
            count++;
            allid++;
        }
        else
        {
            it++;
            allid++;
        }
    }
    return rids;
}

template <int Order>
bool splitIMV(Vector<bool> &ids, Vector<Point> &oldpts, const Curve<2, Order> &crv, const Vector<Real> &dist, Real highBound)
{
    int num = oldpts.size();
    auto polys = crv.getPolys();
    auto knots = crv.getKnots();
    Vector<Point> res(2 * num);
    Vector<bool> resid(2 * num);
    Polynomial<Order, Point> lpoly;

    auto addpt = [&](int &count, Point pt, bool type)
    {
        res[count] = pt;
        resid[count] = type;
        count++;
        if (count >= (int)res.size())
        {
            res.resize(2 * count);
            resid.resize(2 * count);
        }
        return;
    };

    int count = 0;
    addpt(count, oldpts[0], ids[0]);

    for (int i = 0; i < num - 1; i++)
    {
        if (dist[i] <= highBound)
        {
            addpt(count, oldpts[i + 1], ids[i + 1]);
            continue;
        }

        int N = ceil(dist[i] / highBound);
        Real dt = (knots[i + 1] - knots[i]) / N;
        lpoly = polys[i];

        for (int j = 1; j < N; j++)
        {
            addpt(count, lpoly(dt * j), true);
        }
        addpt(count, oldpts[i + 1], ids[i + 1]);
    }
    res.resize(count);
    resid.resize(count);
    oldpts = res;
    ids = resid;
    return count != num;
}

template <int Order>
Vector<unsigned int> MARS2DIMV<Order>::splitLongEdges(const VectorFunction<2> &v, Vector<Point> &pts, const Crv &crv, Real tn, Real dt)
{
    int num = pts.size();
    auto polys = crv.getPolys();
    Vector<Point> oldpts(num);

    for (int i = 0; i < num - 1; i++)
    {
        oldpts[i] = polys[i][0];
    }
    oldpts[num - 1] = polys[0][0];

    Vector<bool> ids(num, false);
    Vector<Real> dist(num - 1);

    while (true)
    {
        for (int i = 0; i < (int)dist.size(); i++)
        {
            dist[i] = norm(pts[i + 1] - pts[i], 2);
        }
        if (!splitIMV(ids, oldpts, crv, dist, (chdLenRange.hi())[0]))
            break;
        pts = oldpts;
        Base::TI->timeStep(v, pts, tn, dt);
        dist.resize(pts.size() - 1);
    }

    Vector<unsigned int> aids(ids.size() - num);
    int count = 0;
    for (int i = 0; i < (int)ids.size(); i++)
    {
        if (ids[i] == true)
        {
            aids[count] = i;
            count++;
        }
    }
    return aids;
}

/*
template <int Order>
Vector<unsigned int> MARS2DIMV<Order>::splitLongEdges(const VectorFunction<2> &v, Vector<Point> &pts, const Crv &crv, Real tn, Real dt)
{
    assert(crv.isClosed(tol));

    int num = pts.size();
    Vector<Point> res(2 * num);
    Vector<unsigned int> ids;
    Vector<Polynomial<Order, Point>> polys = crv.getPolys();
    Polynomial<Order, Point> lpoly;
    Vector<Real> knots = crv.getKnots();
    function<void(Real, Real, Vector<Real> &)> split;

    //lambda: split between two points
    split = [&](Real left, Real right, Vector<Real> &chdt)
    {
        Point lpt = lpoly(left);
        lpt = Base::TI->timeStep(v, lpt, tn, dt);
        Point rpt = lpoly(right);
        rpt = Base::TI->timeStep(v, rpt, tn, dt);
        if (norm(rpt - lpt) > (chdLenRange.hi())[0])
        {
            int N = (int)ceil(norm(rpt - lpt) / (chdLenRange.hi())[0]);
            Real dl = (right - left) / N;
            for (int j = 0; j < N - 1; j++)
            {
                split(left + dl * j, left + dl * (j + 1), chdt);
                chdt.push_back(left + dl * (j + 1));
            }
            split(right - dl, right, chdt);
        }
        return;
    };

    auto addpt = [&](int &count, Point pt)
    {
        res[count] = pt;
        count++;
        if (count >= (int)res.size())
            res.resize(2 * count);
        return;
    };

    unsigned int id = 0;
    int count = 0;
    addpt(count, pts[0]);
    Real dist;
    Point opt;

    for (int i = 0; i < num - 1; i++)
    {
        dist = norm(pts[i + 1] - pts[i]);
        if (dist <= (chdLenRange.hi())[0])
        {
            addpt(count, pts[i + 1]);
            id++;
            continue;
        }

        lpoly = polys[i];
        Vector<Real> chdt; //save parameter t of the added points
        split(0, knots[i + 1] - knots[i], chdt);

        for (int j = 0; j < (int)chdt.size(); j++)
        {
            id++;
            ids.push_back(id);
            opt = lpoly(chdt[j]);
            addpt(count, Base::TI->timeStep(v, opt, tn, dt));
        }
        id++;
        addpt(count, pts[i + 1]);
    }
    res.resize(count);
    pts = res;
    return ids;
}
*/

template <int Order>
void MARS2DIMV<Order>::timeStep(const VectorFunction<2> &v, YS &ys, Real tn, Real dt)
{
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    int id = 1;
    for (auto &crv : vcrv)
    {
        assert(crv.isClosed(tol));

        //pts: now's points
        Vector<Polynomial<Order, Point>> polys = crv.getPolys();
        Vector<Point> pts(polys.size() + 1);
        for (int i = 0; i < (int)polys.size(); i++)
        {
            pts[i] = polys[i][0];
        }
        pts[polys.size()] = polys[0][0];

        //get the points after discrete flow map
        discreteFlowMap(v, pts, tn, dt);

        //split long edges
        Vector<unsigned int> splitids = splitLongEdges(v, pts, crv, tn, dt);

        //remove small edges
        Vector<unsigned int> removeids = removeSmallEdges(pts);

        //
        int num = pts.size();
        Vector<Real> dist(num - 1);

        //use Container<Point> pts to generate the Vector<Point> pts
        //fitCurve only support Vector<Point> version
        for (int i = 0; i < num - 1; i++)
        {
            dist[i] = norm(pts[i + 1] - pts[i], 2);
        }
        crv = fitCurve<Order>(pts, true);

        auto maxp = max_element(dist.begin(), dist.end());
        auto minp = min_element(dist.begin(), dist.end());

        cout << "Curve " << id
             << ":  Add " << splitids.size() << " Points,"
             << " remove " << removeids.size() << " Points."
             << " Max chdlength: " << *maxp
             << " Min chdlength: " << *minp << endl;
        id++;
    }
    ys = YS(SegmentedRealizableSpadjor<Order>(vcrv), tol);
    return;
}

template class MARS2DIMV<4>;