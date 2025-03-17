#include "HighPrecisionNumber.h"
#include <queue>


template <typename Container>
concept HasBeginEnd = requires(Container c) {
  { std::begin(c) } -> std::same_as<typename Container::iterator>;
  { std::end(c) } -> std::same_as<typename Container::iterator>;
};

template <typename T, typename U, typename R, typename Op>
concept BinaryOperation = requires(T lhs, U rhs, R, Op op) {
  { op(lhs, rhs) } -> std::convertible_to<R>;
};

template <HasBeginEnd Container, typename Op = std::plus<void>>
  requires BinaryOperation<typename Container::value_type,
                           typename Container::value_type,
                           typename Container::value_type, Op>
auto highPrecisionAdd(const Container &c, const Op &op = std::plus()) ->
    typename Container::value_type {
  using T = typename Container::value_type;
  if (c.empty()) {
    return T{0};
  };
  auto absCmp = [](const T &lhs, const T &rhs) {
    return fabs(lhs) > fabs(rhs);
  };
  std::priority_queue<T, std::vector<T>, decltype(absCmp)> que;
  for (const auto &e : c) {
    que.push(e);
  }
  while (que.size() > 1) {
    T a = que.top();
    que.pop();
    T b = que.top();
    que.pop();
    que.push(op(a, b));
  }
  return que.top();
}