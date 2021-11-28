#include "MARS2DLIST.h"
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

template <class T>
using List = list<T>;

template <int Order>
void MARS2DLIST<Order>::discreteFlowMap(const VectorFunction<2> &v, List<Point> &pts, Real tn, Real dt)
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
bool removeList(int &count, List<Point> &pts, Real lowBound)
{
    //set distances
    int num = pts.size();
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

template <int Order>
int MARS2DLIST<Order>::removeSmallEdges(List<Point> &pts)
{
    int count = 0;

    while (true)
    {
        if (!removeList(count, pts, (chdLenRange.lo())[0]))
            break;
    }

    return count;
}

template <int Order>
int MARS2DLIST<Order>::splitLongEdges(const VectorFunction<2> &v, List<Point> &pts, const Crv &crv, Real tn, Real dt)
{
    assert(crv.isClosed(1e-15));

    int num = pts.size();
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
    int count = 0;
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
            continue;
        }

        prept = *it;
        Vector<Real> chdt; //save parameter t of the added points
        split(knots[i], knots[i + 1], chdt);

        for (int j = chdt.size() - 1; j >= 0; j--)
        {
            count++;
            Point opt = crv(chdt[j]);
            it = pts.emplace(it, Base::TI->timeStep(v, opt, tn, dt));
        }
        for (int j = 0; j < (int)chdt.size() + 1; j++)
            ++it;
    }
    return count;
}

template <int Order>
void MARS2DLIST<Order>::timeStep(const VectorFunction<2> &v, YS &ys, Vector<List<Point>> &vpts, Real tn, Real dt)
{
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    for (int id = 1; id < (int)vcrv.size() + 1; id++)
    {
        assert(vcrv[id-1].isClosed(1e-15));

        //get the points after discrete flow map
        discreteFlowMap(v, vpts[id-1], tn, dt);

        //split long edges
        int splitids = splitLongEdges(v, vpts[id-1], vcrv[id-1], tn, dt);

        //remove small edges
        int removeids = removeSmallEdges(vpts[id-1]);

        //
        int num = vpts[id-1].size();
        Vector<Real> dist(num - 1);

        //use Container<Point> pts to generate the Vector<Point> pts
        //fitCurve only support Vector<Point> version
        Vector<Point> ptsv(num);
        auto it = vpts[id-1].begin();
        Point prept = *it;
        for (int i = 0; i < num - 1; i++)
        {
            ptsv[i] = prept;
            ++it;
            dist[i] = norm(*it - prept);
            prept = *it;
        }
        ptsv[num - 1] = prept;
        vcrv[id-1] = fitCurve<Order>(ptsv, true);

        auto maxp = max_element(dist.begin(), dist.end());
        auto minp = min_element(dist.begin(), dist.end());

        cout << "Curve " << id
             << ":  Add " << splitids << " Points,"
             << " remove " << removeids << " Points."
             << " Max chdlength: " << *maxp
             << " Min chdlength: " << *minp << endl;
    }
    ys = YS(SegmentedRealizableSpadjor<Order>(vcrv), 1e-15);
    return;
}

template <int Order>
void MARS2DLIST<Order>::trackInterface(const VectorFunction<2> &v, YS &ys, Real StartTime, Real dt, Real EndTime, bool output, string fName, int opstride)
{
    Vector<List<Point>> vpts;
    Vector<Crv> vcrv = ys.getBoundaryCycles();
    for (auto &crv : vcrv)
    {
        Vector<Polynomial<Order, Point>> polys = crv.getPolys();
        List<Point> pts(polys.size() + 1);
        auto it = pts.begin();
        for (int i = 0; i < (int)polys.size(); i++)
        {
            *it = polys[i][0];
            ++it;
        }
        *it = polys[0][0];
        vpts.push_back(pts);
    }

    Real T = StartTime;
    Real t = dt;
    int step = 1;
    if (output == true)
    {
        ofstream of(string(fName + "_Start.dat"), ios_base::binary);
        ys.dump(of);
    }

    while (T < EndTime)
    {
        if (EndTime - T < dt)
        {
            t = EndTime - T;
        }
        cout << "Step: " << step << "     timestep: " << t << endl;

        timeStep(v, ys, vpts, T, t);

        cout << endl;
        T += t;

        if (output == true && T == EndTime)
        {
            ofstream of(string(fName + "_End.dat"), ios_base::binary);
            ys.dump(of);
            break;
        }
        if (output == true && step % opstride == 0)
        {
            ostringstream tmps;
            tmps << step;
            ofstream of(string(fName + "_Step" + tmps.str() + ".dat"), ios_base::binary);
            ys.dump(of);
        }
        step++;
    }
    return;
}

template class MARS2DLIST<4>;