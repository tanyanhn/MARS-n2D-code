#ifndef TESTORIENTEDJORDANCURVE_TY_H
#define TESTORIENTEDJORDANCURVE_TY_H

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class TestOrientedJordanCurve : public CppUnit::TestFixture {
 public:
  void setUp() {}
  void tearDown() {}

  void doTest(int num);
  void test1() { doTest(1); }

  CPPUNIT_TEST_SUITE(TestOrientedJordanCurve);
  CPPUNIT_TEST(test1);
  CPPUNIT_TEST_SUITE_END();
};

#endif  // TESTORIENTEDJORDANCURVE_TY_H