#ifndef TESTORIENTEDJORDANCURVE_TY_H
#define TESTORIENTEDJORDANCURVE_TY_H

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "YinSet/OrientedJordanCurve.h"

class TestOrientedJordanCurve : public CppUnit::TestFixture {
 public:
  void setUp() {}
  void tearDown() {}

  void doTest(int num);
  void test1() { doTest(1); }
  void test2() { doTest(2); }
  void test3() { doTest(3); }

  CPPUNIT_TEST_SUITE(TestOrientedJordanCurve);
  CPPUNIT_TEST(test1);
  CPPUNIT_TEST(test2);
  CPPUNIT_TEST(test3);
  CPPUNIT_TEST_SUITE_END();

 public:
  using rVec = Vec<Real, 2>;
  bool testCircle(const std::string& input, const Real tol,
                  std::string& message);
  bool testRectangle(const std::string& input, const Real tol,
                     std::string& message);
  bool testOrientedJordanCurve(const std::string& input, const Real tol,
                               std::string& message);
  bool verifySpline(const std::vector<Polynomial<4, rVec>>& polys,
                    const std::vector<Real>& knots, const bool periodic,
                    const Real tol);
};

#endif  // TESTORIENTEDJORDANCURVE_TY_H
