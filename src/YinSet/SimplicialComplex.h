#ifndef SIMPLICIALCOMPLEX_TY
#define SIMPLICIALCOMPLEX_TY

#include <algorithm>
#include <cassert>
#include <ostream>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::set;
using std::unordered_map;
using std::unordered_set;
using std::vector;

struct Simplex {
  using Vertex = unsigned long;
  set<Vertex> vertices;

  // constructor
  Simplex() {}

  template <class Containor>
  explicit Simplex(const Containor& v) : vertices(v.begin(), v.end()) {}

  template <typename InputIterator>
  Simplex(InputIterator first, InputIterator last) : vertices(first, last) {}

  // accessor
  int getDimension() const { return vertices.size() - 1; }

  // operator
  bool operator<(const Simplex& rhs) const {
    auto lIt = vertices.begin(), rIt = rhs.vertices.begin(),
         lend = vertices.end(), rend = rhs.vertices.end();
    while (lIt != lend && rIt != rend) {
      if (*lIt != *rIt)
        return *lIt < *rIt;
      ++lIt, ++rIt;
    }
    return rIt != rend;
  }

  bool operator>(const Simplex& rhs) const { return rhs < *this; }

  bool operator==(const Simplex& rhs) const {
    auto lIt = vertices.begin(), rIt = rhs.vertices.begin(),
         lend = vertices.end(), rend = rhs.vertices.end();
    while (lIt != lend && rIt != rend) {
      if (*lIt != *rIt)
        return false;
      ++lIt, ++rIt;
    }
    return true;
  }

  // visualization
  void print(std::ostream& os) const {
    os << "(";
    for (auto i : vertices) {
      os << i << ",";
    }
    os << "), ";
  }
};

// For hash and map
template <>
class std::hash<Simplex> {
  static constexpr std::hash<unsigned int> intHash = {};

 public:
  std::size_t operator()(const Simplex& s) const noexcept {
    size_t res = 0;
    int n = s.getDimension();
    assert(n >= 0 && "Simplex in hash() should be initialed.");
    for (auto i : s.vertices)
      res ^= (intHash(i) << 1);
    return res;
  }
};

template <>
struct std::less<set<Simplex>::iterator> {
  bool operator()(const set<Simplex>::iterator& lhs,
                  const set<Simplex>::iterator& rhs) const {
    return *lhs < *rhs;
  }
};

//======================================================================

class SimplicialComplex {
 protected:
  using SimplexIter = set<Simplex>::iterator;
  using Vertex = Simplex::Vertex;
  vector<set<Simplex>> simplexes;
  unordered_map<Vertex, set<SimplexIter>> mVertex2Simplex;

 public:
  // constructor
  SimplicialComplex(){};

  template <class Containor>
  explicit SimplicialComplex(const Containor& sims);

  template <typename InputIterator>
  SimplicialComplex(InputIterator first, InputIterator last);

  SimplicialComplex(const SimplicialComplex& rhs) { *this = rhs; }

  SimplicialComplex& operator=(const SimplicialComplex& rhs) {
    simplexes.clear();
    mVertex2Simplex.clear();
    for (auto& sims : rhs.simplexes) {
      for (auto& s : sims)
        insert(s);
    }
    return *this;
  }

  // accessor
  const vector<set<Simplex>>& getSimplexes() const { return simplexes; }

  int getDimension() const { return simplexes.size() - 1; }

  // get a vertex's star
  int getStarClosure(Vertex p, SimplicialComplex& closure) const;

  int getLink(Vertex p, unordered_set<Vertex>& res) const;

  // insert a Simplex
  int insert(const Simplex& s);
  int insert(Simplex& s);

  // erase a Simplex
  int erase(const Simplex& s);
  int erase(Simplex& s);

 protected:
  // erase a Simplex, do not consider sub simplex. tool for erase
  template <class Containor>
  int eraseExact(const Containor& containor);

  // find all Simplex appear in every element of sims. tool for erase
  void findShare(
      const vector<unordered_map<Vertex, set<SimplexIter>>::iterator>& sims,
      vector<SimplexIter>& shareSim);

 public:
  // determine if contain a Simplex.
  bool contain(const Simplex& s) const {
    return simplexes[s.getDimension()].find(s) !=
           simplexes[s.getDimension()].end();
  }

  // visualization
  void print(std::ostream& os) const {
    os << "{";
    size_t i = 0;
    for (i = 0; i < simplexes.size(); ++i) {
      for (auto& s : simplexes[i]) {
        s.print(os);
        os << " ";
      }
      os << "\n";
    }
    os << "}" << std::endl;
  }
};

template <class Containor>
inline SimplicialComplex::SimplicialComplex(const Containor& sims) {
  for (auto& s : sims) {
    insert(s);
  }
}

template <typename InputIterator>
inline SimplicialComplex::SimplicialComplex(InputIterator first,
                                            InputIterator last) {
  while (first != last) {
    insert(*first++);
  }
}

inline int SimplicialComplex::getStarClosure(Vertex p,
                                             SimplicialComplex& closure) const {
  auto sims = mVertex2Simplex.find(p);
  if (sims == mVertex2Simplex.end())
    return 0;
  for (auto sIt : sims->second) {
    closure.insert(*sIt);
  }
  return 1;
}

inline int SimplicialComplex::getLink(Vertex p,
                                      unordered_set<Vertex>& res) const {
  auto sims = mVertex2Simplex.find(p);
  if (sims == mVertex2Simplex.end())
    return 0;
  for (auto sIt : sims->second) {
    for (auto vertex : sIt->vertices) {
      res.insert(vertex);
    }
  }
  return 1;
}

inline int SimplicialComplex::insert(const Simplex& s) {
  Simplex copy(s);
  return insert(copy);
}
inline int SimplicialComplex::insert(Simplex& s) {
  int sNSim = s.getDimension();
  if (sNSim == -1)
    return 0;
  if (sNSim >= (int)simplexes.size())
    simplexes.resize(sNSim + 1);
  auto pair = simplexes[sNSim].insert(s);
  auto sIt = pair.first;
  auto b = pair.second;
  if (b == false)
    return 0;

  auto vIt = s.vertices.begin();
  size_t vertex;
  while (vIt != s.vertices.end()) {
    vertex = *vIt;
    mVertex2Simplex[vertex].insert(sIt);
    s.vertices.erase(vIt);
    insert(s);
    vIt = s.vertices.insert(vertex).first;
    ++vIt;
  }
  return 1;
}

inline void SimplicialComplex::findShare(
    const vector<unordered_map<Vertex, set<SimplexIter>>::iterator>& sims,
    vector<SimplexIter>& shareSim) {
  std::less<SimplexIter> cmp;
  size_t i = 0, j = 0, n = sims.size();
  if (n == 1) {
    shareSim.insert(shareSim.end(), sims[0]->second.begin(),
                    sims[0]->second.end());
    return;
  }
  vector<set<SimplexIter>::iterator> arr;
  for (i = 0; i < n; ++i) {
    arr.push_back(sims[i]->second.begin());
  }
  i = 0;
  SimplexIter now;
  while (true) {
    now = *(arr[i]);
    j = i;
    ++i;
    i %= n;
    while (i != j) {
      while (arr[i] != sims[i]->second.end() && cmp(*(arr[i]), now)) {
        ++arr[i];
      }
      if (arr[i] == sims[i]->second.end() || *(arr[i]) != now)
        break;

      ++i;
      i %= n;
    }
    if (i == j) {
      shareSim.push_back(now);
      ++arr[i];
    }

    if (arr[i] == sims[i]->second.end())
      break;
  }
}

inline int SimplicialComplex::erase(const Simplex& s) {
  Simplex copy(s);
  return erase(copy);
}
inline int SimplicialComplex::erase(Simplex& s) {
  int sNSim = s.getDimension();
  auto b = simplexes[sNSim].find(s);
  if (b == simplexes[sNSim].end())
    return 0;

  vector<unordered_map<Vertex, set<SimplexIter>>::iterator> sims;
  vector<SimplexIter> shareSim;
  for (auto vertex : s.vertices) {
    sims.push_back(mVertex2Simplex.find(vertex));
  }
  findShare(sims, shareSim);
  eraseExact(shareSim);
  return 1;
}

template <class Containor>
inline int SimplicialComplex::eraseExact(const Containor& containor) {
  for (auto it : containor) {
    for (auto vertex : it->vertices) {
      auto mit = mVertex2Simplex[vertex].find(it);
      mVertex2Simplex[vertex].erase(mit);
    }
    simplexes[it->getDimension()].erase(it);
  }
  while (!simplexes.empty() && simplexes.back().size() == 0)
    simplexes.pop_back();
  return 1;
}

#endif  // !SIMPLICIALCOMPLEX_TY
