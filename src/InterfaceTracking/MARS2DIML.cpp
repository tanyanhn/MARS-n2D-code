#include "MARS2DIML.h"
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
void MARS2DIML<Order>::discreteFlowMap(const VectorFunction<2> &v, Vector<Point> &pts, Real tn, Real dt)
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
bool removeIML(Vector<unsigned int> &ids, Vector<Point> &pts, Real lowBound)
{
    //set distances
    int num = pts.size();
    Vector<Point> res(num);
    Vector<unsigned int> resid(num);

    //function<void(int&, int&)> addpi;

    //lambda: split between two points
    auto addpi = [&](int &count, int &i) -> void
    {
        res[count] = pts[i];
        resid[count] = ids[i];
        count++;
        i++;
        return;
    };

    Vector<Real> dist(num - 1);
    for (int i = 0; i < num - 1; i++)
    {
        dist[i] = norm(pts[i + 1] - pts[i], 2);
    }

    bool predelete = false; //mark whether the pre-point has been erased
    int count = 0;
    res[count] = pts[0];
    resid[count] = ids[0];
    count++;

    int i = 1;
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
            addpi(count, i);
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
                    addpi(count, i);
                    predelete = false;
                }
            }
            else
            {
                addpi(count, i);
                i++;
                predelete = true;
            }
        }
    }
    if (i == num - 1)
    {
        addpi(count, i);
        res.resize(count);
        pts = res;
        resid.resize(count);
        ids = resid;
        return num != count;
    }
    else
    {
        if (predelete == true || dist[i] >= lowBound)
        {
            addpi(count, i);
            addpi(count, i);
            res.resize(count);
            pts = res;
            resid.resize(count);
            ids = resid;
            return num != count;
        }
        i++;
        addpi(count, i);
        return true;
    }
}

template <int Order>
Vector<unsigned int> MARS2DIML<Order>::removeSmallEdges(Vector<Point> &pts)
{
    int num = pts.size();
    Vector<unsigned int> ids(num);
    for (int i = 0; i < num; i++)
    {
        ids[i] = i;
    }

    while (true)
    {
        if (!removeIML(ids, pts, (chdLenRange.lo())[0]))
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
Vector<unsigned int> MARS2DIML<Order>::splitLongEdges(const VectorFunction<2> &v, Vector<Point> &pts, const Crv &crv, Real tn, Real dt)
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

template <int Order>
void MARS2DIML<Order>::timeStep(const VectorFunction<2> &v, YS &ys, Real tn, Real dt)
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

template class MARS2DIML<4>;