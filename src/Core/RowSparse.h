#ifndef ROWSPARSE_H
#define ROWSPARSE_H

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <map>
#include <ostream>
#include <vector>

#include "Core/Linearization.h"
#include "Core/Tensor.h"

template <class T_RowIndex = int, class T_ColIndex = T_RowIndex>
class RowSparse {
 public:
  RowSparse() = default;

  std::size_t insertRow(const T_RowIndex &r,
                        std::initializer_list<T_ColIndex> cols,
                        std::initializer_list<Real> vals);

  std::size_t insertRow(const T_RowIndex &r, int nz, const T_ColIndex *cols,
                        const Real *vals);

  std::size_t getNumRows() const { return rows.size(); }

  // calculate z = alpha * Ax + y
  template <class T1, class T2, class T3>
  void AXPY(Real alpha, const T1 &x, const T2 &y, T3 &z) const;

  // calculate z = alpha * ax + y, compressed row output
  template <class T1, class T2>
  Tensor<Real, 1> AXPY(Real alpha, const T1 &x, const T2 &y) const;

  // iterators
 public:
  class RowSparseIterator {
    typename std::vector<Real>::const_iterator pVal;
    typename std::vector<T_ColIndex>::const_iterator pCol;
    typename std::vector<int>::const_iterator pNz;
    typename std::vector<T_RowIndex>::const_iterator pRw;

   public:
    friend class RowSparse;
    const Real &getValue(int j) const { return *(pVal + j); }
    const T_ColIndex &getColumn(int j) const { return *(pCol + j); }
    const T_RowIndex &getRow() const { return *pRw; }
    int getNz() const { return *pNz; }

    RowSparseIterator &operator++() {
      pVal += getNz();
      pCol += getNz();
      ++pNz;
      ++pRw;
      return *this;
    }
    void operator++(int) { ++(*this); }
    bool operator!=(const RowSparseIterator &rhs) const {
      return this->pRw != rhs.pRw;
    }
    bool operator==(const RowSparseIterator &rhs) const {
      return this->pRw == rhs.pRw;
    }
    std::ptrdiff_t operator-(const RowSparseIterator &rhs) const {
      return this->pRw - rhs.pRw;
    }
  };

  RowSparseIterator cbegin() const {
    RowSparseIterator it;
    it.pVal = values.cbegin();
    it.pCol = columns.cbegin();
    it.pNz = nonZeros.cbegin();
    it.pRw = rows.cbegin();
    return it;
  }

  RowSparseIterator cend() const {
    RowSparseIterator it;
    it.pRw = rows.cend();
    return it;
  }

  template <class T_Less = std::less<T_RowIndex>>
  RowSparseIterator find_fast(const T_RowIndex &ri, const T_Less &ls) const;

  template <class T_Equal = std::less<T_RowIndex>>
  RowSparseIterator find(const T_RowIndex &ri, const T_Equal &ls) const;

  static void linearIn(const RowSparse<T_RowIndex, T_ColIndex> &a_rs,
                       std::vector<char> *buf);
  static void linearOut(const std::vector<char> &buf, size_t *pos,
                        RowSparse<T_RowIndex, T_ColIndex> *rs);

 protected:
  std::vector<Real> values;
  std::vector<T_ColIndex> columns;
  std::vector<int> nonZeros;
  std::vector<int> rowBegin;
  std::vector<T_RowIndex> rows;

 public:
  template <class RowCompare, class ColCompare>
  void sort(const RowCompare &, const ColCompare &);

  friend std::ostream &operator<<(std::ostream &os,
                                  const RowSparse<T_RowIndex, T_ColIndex> &rs) {
    for (auto rowIt = rs.cbegin(); rowIt != rs.cend(); ++rowIt) {
      os << "Dest: " << rowIt.getRow() << std::endl;
      for (int ct = 0; ct < rowIt.getNz(); ++ct) {
        auto src = rowIt.getColumn(ct);
        os << "Src: " << src << " Val: " << rowIt.getValue(ct) << std::endl;
      }
    }
    return os;
  }
};

//============================================================

template <class T_RowIndex, class T_ColIndex>
inline std::size_t RowSparse<T_RowIndex, T_ColIndex>::insertRow(
    const T_RowIndex &r, int nz, const T_ColIndex *cols, const Real *vals) {
  std::size_t a = rows.size();
  rows.push_back(r);
  rowBegin.push_back(values.size());
  nonZeros.push_back(nz);
  values.insert(values.cend(), vals, vals + nz);
  columns.insert(columns.cend(), cols, cols + nz);
  return a;
}

template <class T_RowIndex, class T_ColIndex>
inline std::size_t RowSparse<T_RowIndex, T_ColIndex>::insertRow(
    const T_RowIndex &r, std::initializer_list<T_ColIndex> cols,
    std::initializer_list<Real> vals) {
  return insertRow(r, cols.size(), cols.begin(), vals.begin());
}

template <class T_RowIndex, class T_ColIndex>
template <class T1, class T2, class T3>
inline void RowSparse<T_RowIndex, T_ColIndex>::AXPY(Real alpha, const T1 &x,
                                                    const T2 &y, T3 &z) const {
  for (auto rowIt = cbegin(); rowIt != cend(); ++rowIt) {
    Real Ax = 0.0;
    for (int k = 0; k < rowIt.getNz(); ++k)
      Ax += rowIt.getValue(k) * x(rowIt.getColumn(k));
    z(rowIt.getRow()) = alpha * Ax + y(rowIt.getRow());
  }
}

template <class T_RowIndex, class T_ColIndex>
template <class T1, class T2>
inline Tensor<Real, 1> RowSparse<T_RowIndex, T_ColIndex>::AXPY(
    Real alpha, const T1 &x, const T2 &y) const {
  Tensor<Real, 1> axpy(getNumRows());
  int i = 0;
  for (auto rowIt = cbegin(); rowIt != cend(); ++rowIt, ++i) {
    Real Ax = 0.0;
    for (int k = 0; k < rowIt.getNz(); ++k)
      Ax += rowIt.getValue(k) * x(rowIt.getColumn(k));
    axpy(i) = alpha * Ax + y(rowIt.getRow());
  }
  return axpy;
}

template <class T_RowIndex, class T_ColIndex>
template <class T_Less>
inline auto RowSparse<T_RowIndex, T_ColIndex>::find_fast(const T_RowIndex &ri,
                                                         const T_Less &ls) const
    -> RowSparseIterator {
  auto it = std::lower_bound(rows.cbegin(), rows.cend(), ri, ls);
  if (it == rows.cend() || ls(ri, *it)) return cend();
  auto i = it - rows.cbegin();
  RowSparseIterator rsit;
  rsit.pRw = it;
  rsit.pNz = nonZeros.cbegin() + i;
  rsit.pCol = columns.cbegin() + rowBegin[i];
  rsit.pVal = values.cbegin() + rowBegin[i];
  return rsit;
}

template <class T_RowIndex, class T_ColIndex>
template <class T_Equal>
inline auto RowSparse<T_RowIndex, T_ColIndex>::find(const T_RowIndex &ri,
                                                    const T_Equal &ls) const
    -> RowSparseIterator {
  auto equal = [&](const T_RowIndex x1) { return ls.compare(x1, ri) == 0; };
  auto it = std::find_if(rows.cbegin(), rows.cend(), equal);
  if (it == rows.cend() || ls.compare(ri, *it) != 0) return cend();
  auto i = it - rows.cbegin();
  RowSparseIterator rsit;
  rsit.pRw = it;
  rsit.pNz = nonZeros.cbegin() + i;
  rsit.pCol = columns.cbegin() + rowBegin[i];
  rsit.pVal = values.cbegin() + rowBegin[i];
  return rsit;
}

template <class T_RowIdx, class T_ColIdx>
inline void RowSparse<T_RowIdx, T_ColIdx>::linearIn(
    const RowSparse<T_RowIdx, T_ColIdx> &a_rs, std::vector<char> *buf) {
  LinearizationHelper::linearIn(a_rs.values, buf);
  LinearizationHelper::linearIn(a_rs.columns, buf);
  LinearizationHelper::linearIn(a_rs.nonZeros, buf);
  LinearizationHelper::linearIn(a_rs.rowBegin, buf);
  LinearizationHelper::linearIn(a_rs.rows, buf);
}

template <class T_RowIdx, class T_ColIdx>
inline void RowSparse<T_RowIdx, T_ColIdx>::linearOut(
    const std::vector<char> &buf, size_t *pos,
    RowSparse<T_RowIdx, T_ColIdx> *rs) {
  LinearizationHelper::linearOut(buf, pos, &(rs->values));
  LinearizationHelper::linearOut(buf, pos, &(rs->columns));
  LinearizationHelper::linearOut(buf, pos, &(rs->nonZeros));
  LinearizationHelper::linearOut(buf, pos, &(rs->rowBegin));
  LinearizationHelper::linearOut(buf, pos, &(rs->rows));
}

template <class T_RowIndex, class T_ColIndex>
template <class RowCompare, class ColCompare>
void RowSparse<T_RowIndex, T_ColIndex>::sort(const RowCompare &rowcompare,
                                             const ColCompare &colcompare) {
  // init
  std::map<T_RowIndex, int, RowCompare> rowIdx;
  for (size_t i = 0; i != rows.size(); ++i) {
    rowIdx[rows[i]] = i;
  }
  std::sort(rows.begin(), rows.end(), rowcompare);
  // reset other data.
  std::vector<Real> tempValues;
  std::vector<T_ColIndex> tempColumns;
  std::vector<int> tempNonZeros;
  std::vector<int> tempRowBegin;
  int pole = 0;
  for (size_t i = 0; i != rows.size(); ++i) {
    // get the row information
    int idx = rowIdx[rows[i]];
    int nz = nonZeros[idx];
    int start = rowBegin[idx];
    tempRowBegin.push_back(pole);
    tempNonZeros.push_back(nz);
    // sort the column
    std::vector<Real> rowvalues(nz);
    std::vector<T_ColIndex> rowcols(nz);
    std::copy(columns.data() + start, columns.data() + start + nz,
              rowcols.data());
    std::map<T_ColIndex, int, ColCompare> colIdx;
    for (int j = 0; j != nz; ++j) {
      colIdx[rowcols[j]] = j;
    }
    std::sort(rowcols.begin(), rowcols.end(), colcompare);
    for (int j = 0; j != nz; ++j) {
      rowvalues[j] = values[start + colIdx[rowcols[j]]];
    }
    for (int j = 0; j < nz; ++j) {
      tempValues.push_back(rowvalues[j]);
      tempColumns.push_back(rowcols[j]);
    }
    pole += nz;
  }
  values.swap(tempValues);
  columns.swap(tempColumns);
  nonZeros.swap(tempNonZeros);
  rowBegin.swap(tempRowBegin);
};

#endif  // ROWSPARSE_H
