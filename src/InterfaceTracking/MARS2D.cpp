#include "MARS2D.h"
#include "Core/Polynomial.h"
#include "Core/Tensor.h"
#include <fstream>
#include <sstream>
#include <list>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <iterator>

Real tol = 1e-15;

using namespace std;

using Point = Vec<Real, 2>;

using Crv = Curve<2, 4>;

template <class T>
using Vector = vector<T>;



template <int Order, template <typename...> class Container>
void MARS2D<Order, Container>::discreteFlowMap(const VectorFunction<2> &v, Container<Point> &pts, Real tn, Real dt)
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
template <template <typename...> class Container>
bool remove(Container<unsigned int> &ids, Container<Vec<Real, 2>> &pts, Real lowBound)
{
    //set distances
    int num = ids.size();
    Vector<Real> dist(num - 1);
    auto pit = pts.begin();
    Point prept = *pit;
    for (int i = 0; i < num - 1; i++)
    {
        ++pit;
        dist[i] = norm(*pit - prept);
        prept = *pit;
    }

    bool predelete = false; //mark whether the pre-point has been erased
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

template <int Order, template <typename...> class Container>
Vector<unsigned int> MARS2D<Order, Container>::removeSmallEdges(Container<Point> &pts)
{
    int num = pts.size();
    Container<unsigned int> ids(num);
    auto it = ids.begin();
    for (int i = 0; i < num; i++)
    {
        *it = i;
        ++it;
    }

    while (true)
    {
        if (!remove(ids, pts, (chdLenRange.lo())[0]))
            break;
    }

    Vector<unsigned int> rids(num - ids.size());
    it = ids.begin();
    unsigned int allid = 0;
    int count = 0;
    while (count < (int)rids.size() - 1)
    {
        if (*it != allid)
        {
            rids[count] = allid;
            count++;
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
Vector<unsigned int> MARS2D<Order, Container>::splitLongEdges(const VectorFunction<2> &v, Container<Point> &pts, const Crv &crv, Real tn, Real dt)
{
    assert(crv.isClosed(tol));

    int num = pts.size();
    Vector<unsigned int> ids;
    Vector<Polynomial<Order, Point>> poly = crv.getPolys();
    Vector<Real> knots = crv.getKnots();
    function<void(Real, Real, Vector<Real> &)> split;

    //lambda: split between two points
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
    unsigned int count = 0;
    Point prept = *it;
    ++it;
    Real dist;

    for (int i = 0; i < num - 1; i++)
    {
        dist = norm(*it - prept);
        if (dist <= (chdLenRange.hi())[0])
        {
            prept = *it;
            ++it;
            count++;
            continue;
        }

        prept = *it;
        Vector<Real> chdt; //save parameter t of the added points
        split(knots[i], knots[i + 1], chdt);

        for (int j = chdt.size() - 1; j >= 0; j--)
        {
            count++;
            ids.push_back(count);
            Point opt = crv(chdt[j]);
            it = pts.emplace(it, Base::TI->timeStep(v, opt, tn, dt));
        }
        count++;
        std::advance(it, (int)chdt.size() + 1);
        //for (int j = 0; j < (int)chdt.size() + 1; j++)
        //    ++it;
    }
    return ids;
}

template <int Order, template <typename...> class Container>
void MARS2D<Order, Container>::timeStep(const VectorFunction<2> &v, YS &ys, Real tn, Real dt)
{
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    int id = 1;
    for (auto &crv : vcrv)
    {
        assert(crv.isClosed(tol));

        //pts: now's points
        Vector<Polynomial<Order, Point>> polys = crv.getPolys();
        Container<Point> pts(polys.size() + 1);
        auto it = pts.begin();
        for (int i = 0; i < (int)polys.size(); i++)
        {
            *it = polys[i][0];
            ++it;
        }
        *it = polys[0][0];

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
        Vector<Point> vpts(num);
        it = pts.begin();
        Point prept = *it;
        for (int i = 0; i < num - 1; i++)
        {
            vpts[i] = prept;
            ++it;
            dist[i] = norm(*it - prept);
            prept = *it;
        }
        vpts[num - 1] = prept;
        crv = fitCurve<Order>(vpts, true);

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

template class MARS2D<4, vector>;
template class MARS2D<4, list>;