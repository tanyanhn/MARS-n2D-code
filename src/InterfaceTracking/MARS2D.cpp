#include "MARS2D.h"
#include <fstream>
#include <sstream>
#include <list>
#include <iomanip>
#include <algorithm>
#include <functional>
//#include <cmath>
#include "Core/Polynomial.h"

Real tol = 1e-15;

using namespace std;

template <int Order, template <typename...> class Container>
void MARS2D<Order, Container>::discreteFlowMap(const IT_VectorFunction<2> &v, Container<Point> &pts, Real tn, Real dt)
{
    Base::TI->timeStep(v, pts, tn, dt);
    return;
}

template <template <typename...> class Container>
bool remove(Container<unsigned int> &ids, Container<Vec<Real, 2>> &pts, Real lowBound)
{
    using Point = Vec<Real, 2>;

    int num = ids.size();
    Vector<Real> dist(num - 1);
    //dist.resize(num - 1);
    auto pit = pts.begin();
    Point prept;
    for (int i = 0; i < num - 1; i++)
    {
        prept = *pit;
        ++pit;
        dist[i] = norm(*pit - prept);
    }

    //mark weather the pre-node is candidate to be erased.
    bool predelete = false;
    auto it = ids.begin();
    pit = pts.begin();
    ++it;
    ++pit;

    int i = 1;
    if (dist[0] < lowBound)
    {
        predelete = true;
        it = ids.erase(it);
        pit = pts.erase(pit);
        i++;
    }

    while (i < num - 2)
    {
        if (dist[i] >= lowBound)
        {
            predelete = false;
            ++it;
            ++pit;
            i++;
        }
        else
        {
            if (dist[i - 1] <= dist[i + 1])
            {
                if (predelete == false)
                {
                    it = ids.erase(it);
                    pit = pts.erase(pit);
                    i++;
                    predelete = true;
                }
                else
                {
                    ++it;
                    ++pit;
                    i++;
                    predelete = false;
                }
            }
            else
            {
                it = ids.erase(++it);
                pit = pts.erase(++pit);
                i += 2;
                predelete = true;
            }
        }
    }
    if (i == num - 1)
    {
        return num != (int)ids.size();
    }
    else
    {
        if (predelete == true || dist[i] >= lowBound)
            return num != (int)ids.size();
        it = ids.erase(it);
        pit = pts.erase(pit);
        return true;
    }
}
/*
template <int Dim, template <typename...> class Container>
bool remove_nondebug(Container<unsigned int> &ids, Container<Vec<Real, Dim>> &pts, Real lowBound)
{
    using Point = Vec<Real, Dim>;

    unsigned int &count = *ids.begin();
    int num = pts.size();
    Vector<Real> dist(num - 1);
    //dist.resize(num - 1);
    auto pit = pts.begin();
    Point prept;
    for (int i = 0; i < num - 1; i++)
    {
        prept = *pit;
        ++pit;
        dist[i] = norm(*pit - prept);
    }

    //mark weather the pre-node is candidate to be erased.
    bool predelete = false;
    pit = pts.begin();
    ++pit;

    int i = 1;
    if (dist[0] < lowBound)
    {
        predelete = true;
        count++;
        pit = pts.erase(pit);
        i++;
    }

    while (i < num - 2)
    {
        if (dist[i] >= lowBound)
        {
            predelete = false;
            ++pit;
            i++;
        }
        else
        {
            if (dist[i - 1] <= dist[i + 1])
            {
                if (predelete == false)
                {
                    count++;
                    pit = pts.erase(pit);
                    i++;
                    predelete = true;
                }
                else
                {
                    ++pit;
                    i++;
                    predelete = false;
                }
            }
            else
            {
                count++;
                pit = pts.erase(++pit);
                i += 2;
                predelete = true;
            }
        }
    }
    if (i == num - 1)
    {
        return num != (int)pts.size();
    }
    else
    {
        if (predelete == true || dist[i] >= lowBound)
            return num != (int)pts.size();
        count++;
        pit = pts.erase(pit);
        return true;
    }
}
*/

template <int Order, template <typename...> class Container>
Vector<unsigned int> MARS2D<Order, Container>::removeSmallEdges(Container<Point> &pts)
{
    int num = pts.size();
    Container<unsigned int> ids;
    for (int i = 0; i < num; i++)
    {
        ids.push_back(i);
    }

    while (true)
    {
        if (!remove(ids, pts, (chdLenRange.lo())[0]))
            break;
    }

    Vector<unsigned int> rids(num - ids.size());
    auto it = ids.begin();
    unsigned int allid = 0;
    int rid = 0;
    while (rid < (int)rids.size() - 1)
    {
        if (*it != allid)
        {
            rids[rid] = allid;
            rid++;
            allid++;
        }
        else
        {
            ++it;
            allid++;
        }
    }
    return rids;
}

template <int Order, template <typename...> class Container>
Vector<unsigned int> MARS2D<Order, Container>::splitLongEdges(const IT_VectorFunction<2> &v, Container<Point> &pts, const Crv &crv, Real tn, Real dt)
{
    assert(crv.isClosed(tol));
    int num = pts.size();
    //Vector<Point> oldpts(num);
    Vector<unsigned int> ids;
    Vector<Polynomial<Order, Point>> poly = crv.getPolys();
    Vector<Real> knots = crv.getKnots();

    /*
    for (int i = 0; i < num - 1; i++)
    {
        oldpts[i] = (poly[i])(0);
    }
    oldpts[num - 1] = (poly[0])(0);
    */

    std::function<void(Real, Real, Vector<Real> &)> split;

    split = [&](Real left, Real right, Vector<Real> &chdt) -> void
    {
        Point lpt = crv(left);
        lpt = Base::TI->timeStep(v, lpt, tn, dt);
        Point rpt = crv(right);
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
        else
        {
            return;
        }
    };

    auto it = pts.begin();
    //auto oit = oldpts.begin();
    //++oit;
    unsigned int count = 0;
    Point prept = *it;
    ++it;

    Real dist;
    for (int i = 0; i < num - 1; i++)
    {
        dist = norm(*it - prept);
        if (dist <= (chdLenRange.hi())[0])
        {
            //++oit;
            prept = *it;
            ++it;
            count++;
            continue;
        }
        Vector<Real> chdt;
        split(knots[i], knots[i + 1], chdt);
        for (int j = 0; j < (int)chdt.size(); j++)
        {
            count++;
            ids.push_back(count);
        }
        count++;
        for (int j = chdt.size() - 1; j >= 0; j--)
        {
            Point opt = crv(chdt[j]);
            //oit = oldpts.emplace(oit, opt);
            it = pts.emplace(it, Base::TI->timeStep(v, opt, tn, dt));
        }
        for (int j = 0; j < (int)chdt.size() + 1; j++)
            ++it;
        //oit += chdt.size() + 1;
    }
    //crv = fitCurve<Order>(oldpts, true);
    return ids;
}

template <int Order, template <typename...> class Container>
void MARS2D<Order, Container>::timeStep(const IT_VectorFunction<2> &v, YS &ys, Real tn, Real dt)
{
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    int id = 1;
    for (auto &crv : vcrv)
    {
        assert(crv.isClosed(tol));

        //pts: now time's points
        Container<Point> pts;
        Vector<Polynomial<Order, Point>> polys = crv.getPolys();
        pts.resize(polys.size() + 1);
        auto it = pts.begin();
        for (int i = 0; i < (int)polys.size(); i++)
        {
            *it = polys[i][0];
            ++it;
        }
        *it = polys[0][0];

        //get the next time's points
        Base::TI->timeStep(v, pts, tn, dt);

        //split long edges, and regenerate the old curve
        Vector<unsigned int> splitids = splitLongEdges(v, pts, crv, tn, dt);

        //remove small edges, update pts and old pts simultaneously
        Vector<unsigned int> removeids = removeSmallEdges(pts);

        int num = pts.size();
        Vector<Point> vpts(num);

        it = pts.begin();
        Point prept = *it;
        Vector<Real> dist(num - 1);
        //dist.resize(num - 1);
        for (int i = 0; i < num - 1; i++)
        {
            vpts[i] = prept;
            ++it;
            dist[i] = norm(*it - prept);
            prept = *it;
        }
        vpts[num - 1] = prept;

        auto maxp = std::max_element(dist.begin(), dist.end());
        auto minp = std::min_element(dist.begin(), dist.end());

        crv = fitCurve<Order>(vpts, true);

        std::cout << "Curve " << id << ":  Add " << splitids.size() << " Points, remove " << removeids.size() << " Points. "
                  << "Max chdlength: " << *maxp << "   Min chdlength: " << *minp << std::endl;
        id++;
    }
    ys = YS(SegmentedRealizableSpadjor<Order>(vcrv), tol);
    return;
}

template class MARS2D<4, vector>;
template class MARS2D<4, list>;