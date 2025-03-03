/**
 * @file Linearization.h
 * @author {JiatuYan} ({2513630371@qq.com})
 * @brief General functions for linearIn and linearOut.
 *
 * @copyright Copyright (c) 2023
 *
 */
#ifndef LINEARIZATION_H_
#define LINEARIZATION_H_
#include <fstream>
#include <iostream>
#include <type_traits>
#include <vector>

namespace LinearizationHelper {

template <class T>
concept need_manually_linearize = !std::is_scalar<T>::value;

template <class T>
inline void linearIn(const T &input, std::vector<char> *buf)
  requires std::is_scalar<T>::value
{
  buf->insert(buf->end(), reinterpret_cast<const char *>(&input),
              reinterpret_cast<const char *>(&input) + sizeof(T));
};

template <class T>
inline void linearIn(const T &input, std::vector<char> *buf)
  requires need_manually_linearize<T>
{
  T::linearIn(input, buf);
};

template <class T>
inline void linearIn(const std::vector<T> &input, std::vector<char> *buf) {
  size_t num_elements = input.size();
  buf->insert(buf->end(), reinterpret_cast<const char *>(&num_elements),
              reinterpret_cast<const char *>(&num_elements) + sizeof(size_t));
  for (auto &element : input) {
    linearIn(element, buf);
  }
};

template <class T>
inline void linearOut(const std::vector<char> &buf, size_t *pos, T *res)
  requires std::is_scalar<T>::value
{
  std::copy(buf.begin() + *pos, buf.begin() + *pos + sizeof(T),
            reinterpret_cast<char *>(res));
  *pos += sizeof(T);
};

template <class T>
inline void linearOut(const std::vector<char> &buf, size_t *pos, T *res)
  requires need_manually_linearize<T>
{
  T::linearOut(buf, pos, res);
};

template <class T>
inline void linearOut(const std::vector<char> &buf, size_t *pos,
                      std::vector<T> *res) {
  size_t num_elements;
  linearOut(buf, pos, &num_elements);
  res->resize(num_elements);
  for (size_t i = 0; i != num_elements; ++i) {
    linearOut(buf, pos, &(res->at(i)));
  }
};

template <class T>
void writeData(const std::string &file_name, const T &data) {
  std::vector<char> buffer;
  LinearizationHelper::linearIn(data, &buffer);
  std::ofstream fout(file_name, std::ios::out | std::ios::binary);
  fout.write(buffer.data(), buffer.size());
  fout.close();
}

template <class T>
void readData(const std::string &file_name, T *data) {
  std::ifstream fin(file_name, std::ios::in | std::ios::binary);
  if (!fin) {
    std::cout << "[Linearization] Read data from file error!" << std::endl;
  }
  fin.seekg(0, fin.end);
  auto size = fin.tellg();
  std::vector<char> buffer(size);
  fin.seekg(0);
  fin.read(buffer.data(), size);
  size_t pos = 0;
  LinearizationHelper::linearOut(buffer, &pos, data);
  fin.close();
}
}  // namespace LinearizationHelper
#endif  // LINEARIZATION_H_