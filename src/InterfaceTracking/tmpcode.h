template <int Dim>
bool remove(Vector<unsigned int> &ids, const Vector<Vec<Real, Dim>> &pts, Real lowBound)
{
    size_t num = ids.size();
    Vector<Real> dist(num - 1);
    //dist.resize(num - 1);
    for (size_t i = 0; i < num - 1; i++)
    {
        dist[i] = norm(pts[ids[i + 1]] - pts[ids[i]]);
    }

    //mark weather the pre-node is candidate to be erased.
    bool predelete = false;
    auto it = ids.begin();
    ++it;

    size_t i = 1;
    if (dist[0] < lowBound)
    {
        predelete = true;
        it = ids.erase(it);
        i++;
    }

    while (i < num - 2)
    {
        if (dist[i] >= lowBound)
        {
            predelete = false;
            ++it;
            i++;
        }
        else
        {
            if (dist[i - 1] <= dist[i + 1])
            {
                if (predelete == false)
                {
                    it = ids.erase(it);
                    i++;
                    predelete = true;
                }
                else
                {
                    ++it;
                    i++;
                    predelete = false;
                }
            }
            else
            {
                it = ids.erase(it + 1);
                i += 2;
                predelete = true;
            }
        }
    }
    if (i == num - 1)
    {
        return num != ids.size();
    }
    else
    {
        if (predelete == true || dist[i] >= lowBound)
            return num != ids.size();
        it = ids.erase(it);
        return true;
    }
}

template <int Dim, int Order>
Vector<unsigned int> MARS<Dim, Order>::removeSmallEdges(Vector<Point> &pts)
{
    size_t Num = pts.size();
    Vector<unsigned int> ids(Num);
    for (size_t i = 0; i < Num; i++)
    {
        ids[i] = i;
    }

    while (true)
    {
        if (!remove(ids, pts, (chdLenRange.lo())[0]))
            break;
    }
    Vector<Point> npts;
    npts.resize(ids.size());
    for (size_t i = 0; i < ids.size(); i++)
    {
        npts[i] = pts[ids[i]];
    }
    pts = npts;
    return ids;
}