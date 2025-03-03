#include "Core/Polynomial.h"
#include "testHeader.h"

Real localTol = 1e-12;

TEST_CASE("Polynomial::translate", "[Poly][Op]") {
  constexpr int Order = 4;
  constexpr int l = -10, r = 20;
  Polynomial<Order, Real> f, g;

  constexpr int tests = 100;
  for (int k = 0; k < tests; ++k) {
    for (int i = 0; i < Order; ++i) {
      f[i] = Generator::randomCreateReal<l, r>();
    }
    Real x = Generator::randomCreateReal<l, r>();
    Real x0 = Generator::randomCreateReal<l, r>();

    Real y = x - x0;
    g = f.translate(x0);
    REQUIRE(std::fabs(g(y) - f(x)) < localTol);

    y = x0 - x;
    g = f.translate(x0, true);
    REQUIRE(std::fabs(g(y) - f(x)) < localTol);

    REQUIRE((-f + f * 1)(x) < localTol);
  }
}
