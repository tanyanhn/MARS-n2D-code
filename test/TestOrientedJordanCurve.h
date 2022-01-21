#ifndef TESTORIENTEDJORDANCURVE_TY_H
#define TESTORIENTEDJORDANCURVE_TY_H

#include "YinSet/OrientedJordanCurve.h"
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class TestOrientedJordanCurve : public CppUnit::TestFixture {
public:
  using rVec = Vec<Real,Dim>;

  void setUp() {}
  void tearDown() {}

  void doTest();
  void test1() { doTest(); }

  void testCircle(const std::string input, const Real tol);
  void testRectangle(const std::string& input, const Real tol);
  void testOrientedJordanCurve(const std::string& input, const Real tol);
  bool TestOrientedJordanCurve::verifySpline(const std::vector<Polynomial<Order, rVec>>& polys,
                                             const std::vector<Real>& knots,
                                             const bool periodic,
                                             const Real tol)
  CPPUNIT_TEST_SUITE(TestOrientedJordanCurve);
  CPPUNIT_TEST(test1);
  CPPUNIT_TEST_SUITE_END();
};

#endif  // TESTORIENTEDJORDANCURVE_TY_H
