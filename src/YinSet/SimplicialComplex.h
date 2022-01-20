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
  set<size_t> vertices;

  Simplex() {}
  template <class Containor>
  explicit Simplex(const Containor& v) : vertices(v.begin(), v.end()) {}
  template <typename InputIterator>
  Simplex(InputIterator first, InputIterator last) : vertices(first, last) {}
  int getNSim() const { return vertices.size() - 1; }
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
  void print(std::ostream& os) const {
    os << "(";
    for (auto i : vertices) {
      os << i << ",";
    }
    os << "), ";
  }
};

template <>
class std::hash<Simplex> {
  static constexpr std::hash<unsigned int> intHash = {};

 public:
  std::size_t operator()(const Simplex& s) const noexcept {
    size_t res = 0;
    int n = s.getNSim();
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
  // constexpr bool operator()(const unordered_set<Simplex>::iterator& lhs,
  //                           const unordered_set<Simplex>::iterator& rhs)
  //                           const {
  //   return *lhs < *rhs;
  // }
};

class SimplicialComplex {
 protected:
  //  public:
  using SimplexIt = set<Simplex>::iterator;
  vector<set<Simplex>> simplexes;
  unordered_map<unsigned int, set<SimplexIt>> mVertex2Simplex;

 public:
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
  const vector<set<Simplex>>& getSimplexes() const { return simplexes; }
  int getNSim() const { return simplexes.size() - 1; }
  int starClosure(unsigned int p, SimplicialComplex& Closure) const;
  int link(unsigned int p, unordered_set<unsigned int>& res) const;
  int insert(const Simplex& s);
  int insert(Simplex& s);
  int erase(const Simplex& s);
  int erase(Simplex& s);
  template <class Containor>
  int eraseExact(const Containor& containor);
  void findShare(
      const vector<unordered_map<unsigned int, set<SimplexIt>>::iterator>& sims,
      vector<SimplexIt>& shareSim);
  bool contain(const Simplex& s) const {
    return simplexes[s.getNSim()].find(s) != simplexes[s.getNSim()].end();
  }
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

inline int SimplicialComplex::starClosure(unsigned int p,
                                          SimplicialComplex& Closure) const {
  auto sims = mVertex2Simplex.find(p);
  if (sims == mVertex2Simplex.end())
    return 0;
  for (auto sIt : sims->second) {
    Closure.insert(*sIt);
  }
  return 1;
}

inline int SimplicialComplex::link(unsigned int p,
                                   unordered_set<unsigned int>& res) const {
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
  int sNSim = s.getNSim();
  if (sNSim == -1)
    return 0;
  if (sNSim >= (int)simplexes.size())
    simplexes.resize(sNSim + 1);
  auto [sIt, b] = simplexes[sNSim].insert(s);
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
    const vector<unordered_map<unsigned int, set<SimplexIt>>::iterator>& sims,
    vector<SimplexIt>& shareSim) {
  std::less<SimplexIt> cmp;
  size_t i = 0, j = 0, n = sims.size();
  if (n == 1) {
    shareSim.insert(shareSim.end(), sims[0]->second.begin(),
                    sims[0]->second.end());
    return;
  }
  vector<set<SimplexIt>::iterator> arr;
  for (i = 0; i < n; ++i) {
    arr.push_back(sims[i]->second.begin());
  }
  i = 0;
  SimplexIt now;
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
  int sNSim = s.getNSim();
  auto b = simplexes[sNSim].find(s);
  if (b == simplexes[sNSim].end())
    return 0;

  vector<unordered_map<unsigned int, set<SimplexIt>>::iterator> sims;
  vector<SimplexIt> shareSim;
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
    simplexes[it->getNSim()].erase(it);
  }
  while (!simplexes.empty() && simplexes.back().size() == 0)
    simplexes.pop_back();
  return 1;
}

#endif  // !SIMPLICIALCOMPLEX_TY
