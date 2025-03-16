#ifndef _COMPUTEERROR_H_
#define _COMPUTEERROR_H_

#include <cmath>
#include <cstddef>
#include <numeric>

#include "InterfaceTracking/XorArea.h"
#include "Marsn2D/InterfaceGraph.h"
#include "Recorder/Recorder.h"
#include "YinSet/YinSet.h"

using Crv = Curve<2, 4>;
using Point = Vec<Real, 2>;

template <class T>
using Vector = std::vector<T>;

inline Vector<Real> circleerror(const Vector<Crv> &crvs, const Point &center,
                                Real radio) {
  size_t n = crvs.size();
  Vector<Real> result(2 * n - 1);
  for (int i = 0; i < n; i++) {
    auto polys = crvs[i].getPolys();
    auto knots = crvs[i].getKnots();
    size_t cn = polys.size();
    Vector<Real> dis(cn);
    for (int j = 0; j < cn; j++) {
      auto pt = crvs[i]((knots[j] + knots[j + 1]) / 2);
      dis[j] = abs(norm(pt - center, 2) - radio);
    }

    auto sum = accumulate(dis.begin(), dis.end(), 0.0);
    result[2UL * i] = sum / (cn + 1);

    /*
    std::cout << std::endl;
    auto maxp = max_element(dis.begin(), dis.end());
    result[2 * i] = *maxp;
    */
  }
  for (int i = 0; i < n - 1; i++) {
    result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
  }
  return result;
}

Vector<Real> squarerror(const Vector<YinSet<2, 4>> &vys, const Point &center,
                        Real half) {
  int n = vys.size();
  Vector<Real> result(2 * n - 1);
  for (int i = 0; i < n; i++) {
    auto crv = vys[i].getBoundaryCycles()[0];
    SimplicialComplex<typename YinSet<2, 4>::PointIndex> kinks =
        vys[i].getKinks();
    Vector<unsigned int> brk;
    if (!kinks.getSimplexes().empty()) {
      for (const auto &simplex : kinks.getSimplexes()[0]) {
        typename YinSet<2, 4>::Vertex id = *simplex.vertices.begin();
        brk.push_back(id.second);
      }
    }
    assert(brk[0] == 0 && brk.size() == 4);
    const auto &polys = crv.getPolys();
    const auto &knots = crv.getKnots();
    size_t cn = polys.size();
    Vector<Real> dis(cn);
    for (unsigned int j = brk[0]; j < brk[1]; j++) {
      auto pt = crv(knots[j]);
      dis[j] = abs(pt[1] - center[1] + half);
    }
    for (unsigned int j = brk[1]; j < brk[2]; j++) {
      auto pt = crv(knots[j]);
      dis[j] = abs(pt[0] - center[0] - half);
    }
    for (unsigned int j = brk[2]; j < brk[3]; j++) {
      auto pt = crv(knots[j]);
      dis[j] = abs(pt[1] - center[1] - half);
    }
    for (unsigned int j = brk[3]; j < (unsigned int)cn; j++) {
      auto pt = crv(knots[j]);
      dis[j] = abs(pt[0] - center[0] + half);
    }
    /*
    auto sum = accumulate(dis.begin(), dis.end(), 0.0);
    result[2 * i] = sum / (cn + 1);
    */
    auto maxp = max_element(dis.begin(), dis.end());
    result[2 * i] = *maxp;
  }
  for (int i = 0; i < n - 1; i++) {
    result[2 * i + 1] = log(result[2UL * i] / result[2 * i + 2]) / log(2);
  }
  return result;
}

inline Vector<Real> exactError(const Vector<Crv> &crvs, const Crv &exact,
                               Real tol) {
  int n = crvs.size();
  Vector<Real> result(2 * n - 1);
  for (int i = 0; i < n; i++) {
    result[2 * i] = xorArea(exact, crvs[i], tol);
    // result[2 * i] = xorArea(crvs[i], exact, tol);
  }
  for (int i = 0; i < n - 1; i++) {
    result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
  }
  return result;
}

inline Vector<Real> richardsonError(const Vector<Crv> &crvs, Real tol) {
  int n = crvs.size();
  Vector<Real> result(2 * n - 3);
  for (int i = 0; i < n - 1; i++) {
    // result[2 * i] = xorArea(crvs[i + 1], crvs[i], tol);
    result[2 * i] = xorArea(crvs[i], crvs[i + 1], tol);
  }
  for (int i = 0; i < n - 2; i++) {
    result[2 * i + 1] = log(result[2 * i] / result[2 * i + 2]) / log(2);
  }
  return result;
}
inline void dumpError(const Vector<Vector<Real>> &errors,
                      const std::string &filename) {
  std::ofstream of(filename);
  auto N1 = errors.size();
  auto N2 = errors[0].size();
  of.write((char *)&N1, sizeof(int));
  of.write((char *)&N2, sizeof(int));
  for (int i = 0; i < N1; i++) {
    for (int j = N2 - 1; j >= 0; j--) {
      of.write((char *)&errors[i][j], sizeof(Real));
    }
  }
}

// return: each yinset - 2(L^2 and L^\infty) - each grid
template <int Order>
#ifdef OPTNONE
__attribute__((optnone))
#endif  // OPTNONE
Vector<Vector<Vector<Real>>>
cutCellError(const vector<Marsn2D::approxInterfaceGraph<Order>> &lhss,
             const Marsn2D::approxInterfaceGraph<Order> &rhs, auto &boxs,
             auto &range) {
  const int num = lhss.size();  // number of grid layers
  int numYinsets = 0;           // number of yinsets in each layer
  // compute volumes in a box domain, return Vector<Vector<Real>>
  auto calculateVolume = [range](const auto &box, const auto &yinset) {
    Vec<Real, 2> h;
    h[0] = range.hi()[0] - range.lo()[0];
    h[1] = range.hi()[1] - range.lo()[1];
    h = h / (box.size());
    Real fullCell = h[0] * h[1];
    auto N = box.size();
    Vector<Vector<Real>> volume(N[0], Vector<Real>(N[1]));
    auto [rhsRes, rhsBoundary, rhsTags] = yinset.cutCell(box, range, false);
    // std::ofstream of(getExportDir() + "localYinset" + ".dat",
    // std::ios_base::binary); dumpVecYinSet<Order>(rhsRes, of);
    loop_box_2(box, i0, i1) {
      if (rhsTags[i0][i1] == 1 || rhsTags[i0][i1] == 2) {
        volume[i0][i1] = fullCell;
      } else if (rhsTags[i0][i1] == -1 || rhsTags[i0][i1] == -2) {
        volume[i0][i1] = 0;
      } else if (rhsTags[i0][i1] == 0) {
        volume[i0][i1] = 0;
        if (rhsRes[i0][i1]) {
          for (const auto &crv : rhsRes[i0][i1]->getBoundaryCycles())
            volume[i0][i1] += area(crv);
        }
      } else {
        throw std::runtime_error("Invalid tag");
      }
    }
    return volume;
  };
  auto coarseVolumes = [](Vector<Vector<Vector<Real>>> volumes) {
    const size_t Nx = volumes[0].size();
    const size_t Ny = volumes[0][0].size();
    Vector<Vector<Vector<Real>>> ret(
        volumes.size(), Vector<Vector<Real>>(Nx / 2, Vector<Real>(Ny / 2, 0)));

    for (size_t k = 0; k < volumes.size(); ++k) {
      for (size_t i = 0; i < Nx / 2; i++) {
        for (size_t j = 0; j < Ny / 2; j++) {
          ret[k][i][j] =
              volumes[k][2 * i][2 * j] + volumes[k][2 * i + 1][2 * j] +
              volumes[k][2 * i][2 * j + 1] + volumes[k][2 * i + 1][2 * j + 1];
        }
      }
    }
    return ret;
  };

  // compute rhs
  auto &box = boxs.back();
  const auto &rhsYinsets = rhs.approxYinSet();
  Vector<Vector<Vector<Real>>> rhsVolumes;  // each yinset - each cell (2D)
  numYinsets = rhsYinsets.size();
  for (size_t i = 0; i < rhsYinsets.size(); i++) {
    rhsVolumes.emplace_back(calculateVolume(box, rhsYinsets[i]));
  }

  // compute lhss
  Vector<Vector<Vector<Real>>> ret(
      numYinsets, Vector<Vector<Real>>(2, Vector<Real>(2 * num - 1, 0)));
  for (size_t g = 0; g < num; g++) {
    auto grid = num - g - 1;
    Real h = 1.0 / (1 << grid);
    const auto &lhsYinsets = lhss[grid].approxYinSet();
    for (size_t i = 0; i < numYinsets; ++i) {
      auto &L1 = ret[i][0][2 * grid];
      auto &LInf = ret[i][1][2 * grid];
      auto localVolumes = calculateVolume(boxs[grid], lhsYinsets[i]);
      // if (i == 2) {
      //   dumpError(rhsVolumes[i], getExportDir() + "localVolumes" + ".dat");
      //   dumpError(localVolumes,
      //             getExportDir() + "localVolumes" + std::to_string(i) +
      //             ".dat");
      // }
      loop_box_2(boxs[grid], i0, i1) {
        L1 += std::fabs(localVolumes[i0][i1] - rhsVolumes[i][i0][i1]);
        LInf = std::max(LInf,
                        abs(localVolumes[i0][i1] - rhsVolumes[i][i0][i1]) / h);
        // if (LInf > 1e-6) breakpoint();
      }
    }
    rhsVolumes = coarseVolumes(rhsVolumes);
  }
  for (size_t i = 0; i < num - 1; i++) {
    for (size_t j = 0; j < numYinsets; ++j) {
      ret[j][0][2 * i + 1] =
          log(ret[j][0][2 * i] / ret[j][0][2 * i + 2]) / log(2);
      ret[j][1][2 * i + 1] =
          log(ret[j][1][2 * i] / ret[j][1][2 * i + 2]) / log(2);
    }
  }
  return ret;
}

template <typename Range>
std::string join(const Range &range) {
  std::ostringstream oss;
  bool first = true;
  int count = 0;
  for (const auto &elem : range) {
    if (!first) oss << ' ';
    if (++count % 2)
      oss << std::format("{:.2e} ", elem);
    else
      oss << std::format("{:.2f} ", elem);
    first = false;
  }
  return oss.str();
}

inline void printCellError(const Vector<Vector<Vector<Real>>> &errors) {
  int count = 0;
  for (const auto &i : errors) {
    std::cout << std::format("YinSet-{} L1    : {}\n         Linfty: {}\n\n",
                             count++, join(i[0]), join(i[1]));
  }
}

#endif