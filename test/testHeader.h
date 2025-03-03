// Only include by test file.
#include "Core/Config.h"
#include "Core/dirConfig.h"
#include "catch_amalgamated.hpp"

using namespace std;
using namespace Catch;

struct Generator {
  template <int l, int r>
  [[nodiscard]] static auto randomCreateReal() -> Real;
};

const Generator generator;

template <int l, int r>
auto Generator::randomCreateReal() -> Real {
  STATIC_CHECK(l < r);
  return GENERATE(take(1, random(l, r)));
}