#ifndef TESTYINSETSBOOLEANOPS_H
#define TESTYINSETSBOOLEANOPS_H

#include <string>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "Core/Vec.h"

class TestYinSetsBooleanOps : public CppUnit::TestFixture
{
public:
  using string = std::string;
  using rVec = Vec<Real,2>;
  
  void setUp() { }
  void tearDown() { }

  void doTest(const string &input1, const string &input2, Real tol,
              const rVec &offset, const string &answer);

  void test1() {
    doTest("data/spadjor-3.input.dat","data/spadjor-4.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-3-spadjor-4.answer.dat");
  }
  void test2() {
    doTest("data/spadjor-3.input.dat","data/spadjor-5.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-3-spadjor-5.answer.dat");
  }
  void test3() {
    doTest("data/spadjor-3.input.dat","data/spadjor-6.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-3-spadjor-6.answer.dat");
  }
  void test4() {
    doTest("data/spadjor-3.input.dat","data/spadjor-7.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-3-spadjor-7.answer.dat");
  }
  void test5() {
    doTest("data/spadjor-3.input.dat","data/spadjor-10.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-3-spadjor-10.answer.dat");
  }
  void test6() {
    doTest("data/spadjor-8.input.dat","data/spadjor-9.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-8-spadjor-9.answer.dat");
  }
  void test7() {
    doTest("data/spadjor-11.input.dat","data/spadjor-12.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-11-spadjor-12.answer.dat");
  }
  void test8() {
    doTest("data/spadjor-11.input.dat","data/spadjor-13.input.dat",
	   1e-10, rVec{0,0},
	   "data/spadjor-11-spadjor-13.answer.dat");
  }

  void test101() {
    doTest("data1/Panda.dat","data/mickey.input.dat",
	   1e-10, rVec{-1,0},
	   "data/panda-mickey--1.answer.dat");
  }
  void test102() {
    doTest("data1/Panda.dat","data/mickey.input.dat",
	   1e-10, rVec{0,0},
	   "data/panda-mickey-0.answer.dat");
  }
  void test103() {
    doTest("data1/Panda.dat","data/mickey.input.dat",
	   1e-10, rVec{1,0},
	   "data/panda-mickey-1.answer.dat");
  }

  CPPUNIT_TEST_SUITE(TestYinSetsBooleanOps);
  CPPUNIT_TEST(test1);
  CPPUNIT_TEST(test2);
  CPPUNIT_TEST(test3);
  CPPUNIT_TEST(test4);
  CPPUNIT_TEST(test5);
  CPPUNIT_TEST(test6);
  CPPUNIT_TEST(test7);
  CPPUNIT_TEST(test8);
  CPPUNIT_TEST(test101);
  CPPUNIT_TEST(test102);
  CPPUNIT_TEST(test103);
  CPPUNIT_TEST_SUITE_END();
  
};

#endif // TESTYINSETSBOOLEANOPS_H
