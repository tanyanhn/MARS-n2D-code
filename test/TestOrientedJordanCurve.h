#ifndef TESTORIENTEDJORDANCURVE_TY_H
#define TESTORIENTEDJORDANCURVE_TY_H

#include "YinSet/OrientedJordanCurve.h"
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class TestOrientedJordanCurve : public CppUnit::TestFixture {
public:

  void setUp() {}
  void tearDown() {}

  void doTest();
  void test1() { doTest(); }

  CPPUNIT_TEST_SUITE(TestOrientedJordanCurve);
  CPPUNIT_TEST(test1);
  CPPUNIT_TEST_SUITE_END();
public:
  using rVec = Vec<Real,2>;
  bool testCircle(const std::string& input, const Real tol, std::string& message);
  bool testRectangle(const std::string& input, const Real tol, std::string& message);
  bool testOrientedJordanCurve(const std::string& input, const Real tol, std::string& message);
  bool verifySpline(const std::vector<Polynomial<4, rVec>>& polys,
                    const std::vector<Real>& knots,
                    const bool periodic,
                    const Real tol);
};

#endif  // TESTORIENTEDJORDANCURVE_TY_H
