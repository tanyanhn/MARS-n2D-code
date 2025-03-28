#ifndef TENSOR_H
#define TENSOR_H

#include <iomanip>

#include "Box.h"
#include "Recorder/Recorder.h"
#include "dirConfig.h"

template <class T, int Dim>
class Tensor;
template <class T, int Dim>
class TensorSlice;
template <class, class>
struct TensorUnaryOp;
template <class, class, class>
struct TensorBinaryOp;

//==================================================

template <class TExpr>
class TensorExpr {
 public:
  // i=0 points to the first meaningful element, return non-reference
  template <class... Ts>
  auto at(int i, Ts... args) const {
    return static_cast<const TExpr&>(*this).at(i, args...);
  }
  //
  auto box() const { return static_cast<const TExpr&>(*this).box(); }
};

//==================================================

template <class T, int Dim>
class Tensor : public TensorExpr<const Tensor<T, Dim>&> {
 public:
  using iVec = Vec<int, Dim>;

  Tensor() : aData(nullptr), /*zData(nullptr), */ owned(false) {}

  Tensor(const Box<Dim>& _box);

  Tensor(const Box<Dim - 1>& _box, int ncomp)
      : Tensor(enlarge(_box, Box<1>(0, ncomp - 1))) {
    //    static_assert(Dim>=2, "invalid dimension.");
  }

  Tensor(const iVec& sz) : Tensor(Box<Dim>(0, sz - 1)) {}

 public:
  // ========================================= copy ctor begin
  // copy ctor is for aliasing
  Tensor(const Tensor<T, Dim>& rhs);

  Tensor(const Tensor<T, Dim + 1>& rhs, int comp);

  // copy construct from expression
  template <class TExpr>
  Tensor(const TensorExpr<TExpr>& expr);

  // ========================================= copy assignment begin
  // copy assignment is for deep copy
  Tensor<T, Dim>& operator=(const Tensor<T, Dim>& rhs);

  // calculate results from expression
  template <class TExpr>
  Tensor<T, Dim>& operator=(const TensorExpr<TExpr>& expr);

  // ========================================= move ctor and move assignment
  // begin move ctor is for stealing
  Tensor(Tensor<T, Dim>&& rhs) noexcept;
  // move assignment is for stealing
  Tensor<T, Dim>& operator=(Tensor<T, Dim>&& rhs) noexcept;

  // ========================================= other operations begin
  // force copy
  Tensor<T, Dim> copy() const;

  void resize(const Box<Dim>& _box);
  void resize(const Box<Dim - 1>& _box, int ncomp);
  void resize(const iVec& sz);
  void clear();
  ~Tensor() { clear(); }

  template <class Q, int D>
  friend void swap(Tensor<Q, D>& lhs, Tensor<Q, D>& rhs);

  void dump(std::ostream& os) const;
  void dump(const char* file) const;

  void load(const char* file_name);

  //========================================
  // initialization
 public:
  Tensor<T, Dim>& operator=(std::initializer_list<T> initlist);
  Tensor<T, Dim>& operator=(const T& a);

  //========================================
  // slicing
 public:
  TensorSlice<T, Dim> slice(const Box<Dim>& pbox) const;

  // Note : This is a co-dimension 1 reduction.
  // Use this syntax repeatedly for reduction of higher co-dimensions.
  TensorSlice<T, Dim - 1> slice(int D, int a) const;

  //========================================
  // getters
 public:
  const T* __restrict__ data() const { return aData; }
  T* __restrict__ data() { return aData; }

  const Box<Dim>& box() const { return bx; }

  const iVec& getStride() const { return stride; }

  iVec size() const { return bx.size(); }

  int volume() const { return bx.volume(); }

  void view(std::ostream& os = std::cout, int precision = 3) const {
    std::streamsize origin = os.precision();
    os << std::setprecision(precision) << std::scientific;
    os << *this << std::endl;
    os.precision(origin);
    os << std::defaultfloat;
  }

  void view(const char* file, int precision = 2) const {
    std::ofstream ofs(file);
    ofs << std::setprecision(precision);
    ofs << *this << std::endl;
    ofs.close();
  }

  //========================================
  // naive accessor
 public:
  const T& operator()(const Vec<int, Dim>& idx) const {
    return aData[dot(idx - bx.lo(), stride)];
  }
  T& operator()(const Vec<int, Dim>& idx) {
    return aData[dot(idx - bx.lo(), stride)];
  }
  // one dimension
  const T& operator()(int i0) const { return aData[i0 - bx.lo()[0]]; }
  T& operator()(int i0) { return aData[i0 - bx.lo()[0]]; }
  // two dimensions
  const T& operator()(int i0, int i1) const {
    return aData[(i0 - bx.lo()[0]) + (i1 - bx.lo()[1]) * stride[1]];
  }
  T& operator()(int i0, int i1) {
    return aData[(i0 - bx.lo()[0]) + (i1 - bx.lo()[1]) * stride[1]];
  }
  // three dimensions
  const T& operator()(int i0, int i1, int i2) const {
    return aData[(i0 - bx.lo()[0]) + (i1 - bx.lo()[1]) * stride[1] +
                 (i2 - bx.lo()[2]) * stride[2]];
  }
  T& operator()(int i0, int i1, int i2) {
    return aData[(i0 - bx.lo()[0]) + (i1 - bx.lo()[1]) * stride[1] +
                 (i2 - bx.lo()[2]) * stride[2]];
  }
  // four dimensions
  const T& operator()(int i0, int i1, int i2, int i3) const {
    return aData[(i0 - bx.lo()[0]) + (i1 - bx.lo()[1]) * stride[1] +
                 (i2 - bx.lo()[2]) * stride[2] + (i3 - bx.lo()[3]) * stride[3]];
  }
  T& operator()(int i0, int i1, int i2, int i3) {
    return aData[(i0 - bx.lo()[0]) + (i1 - bx.lo()[1]) * stride[1] +
                 (i2 - bx.lo()[2]) * stride[2] + (i3 - bx.lo()[3]) * stride[3]];
  }

  // support for template expressions
 protected:
  template <class TExpr>
  friend class TensorExpr;
  template <class, class>
  friend struct TensorUnaryOp;
  template <class, class, class>
  friend struct TensorBinaryOp;
  // let i=0 points to the first meaningful element
  const T& at(int i) const { return aData[i]; }
  T& at(int i) { return aData[i]; }
  const T& at(int i, int j) const { return aData[i + j * stride[1]]; }
  T& at(int i, int j) { return aData[i + j * stride[1]]; }
  const T& at(int i, int j, int k) const {
    return aData[i + j * stride[1] + k * stride[2]];
  }
  T& at(int i, int j, int k) {
    return aData[i + j * stride[1] + k * stride[2]];
  }
  const T& at(int i, int j, int k, int l) const {
    return aData[i + j * stride[1] + k * stride[2] + l * stride[3]];
  }
  T& at(int i, int j, int k, int l) {
    return aData[i + j * stride[1] + k * stride[2] + l * stride[3]];
  }

  static iVec strideFromSize(const iVec& sz);

 protected:
  T* __restrict__ aData;  // point to bx.lo()
                          //  T         *zData; // point to (0,...,0)
  bool owned;
  iVec stride;
  Box<Dim> bx;
};

//================================================

template <class T, int Dim>
inline Vec<int, Dim> Tensor<T, Dim>::strideFromSize(const Vec<int, Dim>& sz) {
  Vec<int, Dim> stride;
  stride[0] = 1;
  for (int d = 1; d < Dim; d++) stride[d] = stride[d - 1] * sz[d - 1];
  return stride;
}

template <class T, int Dim>
inline Tensor<T, Dim>::Tensor(const Box<Dim>& _box) : bx(_box) {
  stride = strideFromSize(bx.size());
  aData = new T[volume()];
  //  zData = aData - dot(bx.lo(), stride);
  owned = true;
}

template <class T, int Dim>
inline Tensor<T, Dim>::Tensor(const Tensor<T, Dim>& rhs)
    : owned(false), stride(rhs.stride), bx(rhs.bx) {
  aData = const_cast<T*>(rhs.aData);
  //  zData = const_cast<T*>(rhs.zData);
}

template <class T, int Dim>
inline Tensor<T, Dim>::Tensor(const Tensor<T, Dim + 1>& rhs, int comp) {
  const auto& st = rhs.getStride();
  stride = reduce(st);
  aData = const_cast<T*>(rhs.data() + comp * st[Dim]);
  //  zData = const_cast<T*>(rhs.zdata() + comp * st[Dim]);
  bx = reduce(rhs.box());
  owned = false;
}

template <class T, int Dim>
inline Tensor<T, Dim>& Tensor<T, Dim>::operator=(const Tensor<T, Dim>& rhs) {
  assert(volume() == rhs.volume());
  int vol = volume();
#ifdef USE_TENSOR_OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (int i = 0; i < vol; ++i) aData[i] = rhs.aData[i];
  return *this;
}

template <class T, int Dim>
inline Tensor<T, Dim>::Tensor(Tensor<T, Dim>&& rhs) noexcept : Tensor() {
  swap(*this, rhs);
}

template <class T, int Dim>
inline Tensor<T, Dim>& Tensor<T, Dim>::operator=(
    Tensor<T, Dim>&& rhs) noexcept {
  swap(*this, rhs);
  return *this;
}

//==================================================

template <class T, int Dim>
inline Tensor<T, Dim> Tensor<T, Dim>::copy() const {
  Tensor<T, Dim> res(box());
  res = *this;
  return res;
}

template <class T, int Dim>
inline void Tensor<T, Dim>::resize(const Box<Dim>& _box) {
  clear();
  new (this) Tensor(_box);
}

template <class T, int Dim>
inline void Tensor<T, Dim>::resize(const Box<Dim - 1>& _box, int ncomp) {
  clear();
  new (this) Tensor(_box, ncomp);
}

template <class T, int Dim>
inline void Tensor<T, Dim>::resize(const iVec& sz) {
  resize(Box<Dim>(0, sz - 1));
}

template <class T, int Dim>
inline void Tensor<T, Dim>::clear() {
  if (owned) delete[] aData;
  aData = /*zData = */ nullptr;
  owned = false;
}

template <class T, int Dim>
inline void swap(Tensor<T, Dim>& lhs, Tensor<T, Dim>& rhs) {
  std::swap(lhs.aData, rhs.aData);
  //  std::swap(lhs.zData, rhs.zData);
  std::swap(lhs.owned, rhs.owned);
  std::swap(lhs.stride, rhs.stride);
  std::swap(lhs.bx, rhs.bx);
}

//================================================
// initialization

template <class T, int Dim>
inline Tensor<T, Dim>& Tensor<T, Dim>::operator=(
    std::initializer_list<T> initlist) {
  int i = 0;
  for (auto j = initlist.begin(); j != initlist.end(); ++j, ++i) aData[i] = *j;
  return *this;
}

template <class T, int Dim>
inline Tensor<T, Dim>& Tensor<T, Dim>::operator=(const T& a) {
  int vol = volume();
#ifdef USE_TENSOR_OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
  for (int i = 0; i < vol; ++i) aData[i] = a;
  return *this;
}

template <class T, int Dim>
inline void Tensor<T, Dim>::dump(std::ostream& os) const {
  int z = Dim;
  os.write((char*)&z, sizeof(int));
  const auto sz = size();
  os.write((char*)sz.data(), sizeof(int) * Dim);
  os.write((char*)aData, sizeof(T) * volume());
}

template <class T, int Dim>
void Tensor<T, Dim>::dump(const char* name) const {
  std::string file = getExportDir() + "/" + name;
  std::ofstream os(file, std::ios::binary);
  int z = Dim;
  os.write((char*)&z, sizeof(int));
  const auto sz = size();
  os.write((char*)sz.data(), sizeof(int) * Dim);
  os.write((char*)aData, sizeof(T) * volume());
}

template <class T, int Dim>
inline void Tensor<T, Dim>::load(const char* file_name) {
  std::string file = std::string(ROOT_DIR) + file_name;
  std::ifstream fin(file, std::ios::binary);
  if (!fin.is_open()) {
    std::cout << "Cannot open file: " << file << std::endl;
  } else {
    int d;
    Vec<int, Dim> sz;
    fin.read((char*)&d, sizeof(int));
    fin.read((char*)sz.data(), sizeof(int) * Dim);
    resize(sz);
    fin.read((char*)aData, sizeof(T) * volume());
  }
  fin.close();
}

//================================================
// formatting

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Tensor<T, 1>& t) {
  os << "(";
  loop_box_1(t.box(), i) os << t(i) << " ";
  os << ")";
  return os;
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Tensor<T, 2>& t) {
  os << "(\n";
  const auto& box = t.box();
  for (int i = box.lo()[0]; i <= box.hi()[0]; ++i) {
    for (int j = box.lo()[1]; j <= box.hi()[1]; ++j)
      os << std::setprecision(2) << std::scientific << t(i, j) << " ";
    os << "\n";
  }
  os << ")";
  return os;
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Tensor<T, 3>& t) {
  os << "{\n";
  auto comps = getComp(t.box(), 2);
  loop_box_1(comps, k) {
    Tensor<T, 2> a(t, k);
    os << a;
  }
  os << "}";
  return os;
}

template <class T, int Dim>
void dumpTensor(const char* root, const char* name, const Tensor<T, Dim>& t) {
  // root must begin with '/' and end with '/', relative to ROOT_DIR
  std::ofstream fout;
  std::string file_name = std::string(root) + "/" + name;
  fout.open(file_name);
  if (!fout.is_open()) {
    std::cout << "Cannot open file: " << file_name << std::endl;
  } else {
    t.dump(fout);
  }
  fout.close();
}

template <class T, int Dim>
void dumpDebugTensor(const char* name, const Tensor<T, Dim>& t) {
  dumpTensor<T, Dim>(getExportDir().c_str(), name, t);
}

#endif  // TENSOR_H
