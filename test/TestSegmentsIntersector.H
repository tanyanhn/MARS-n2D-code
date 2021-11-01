#ifndef TESTSEGMENTSINTERSECTOR_H
#define TESTSEGMENTSINTERSECTOR_H

#include "YinSet/SegmentsIntersector.h"
#include <string>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

using std::string;
using std::vector;

class TestSegmentsIntersector : public CppUnit::TestFixture
{
public:
  void setUp() { }
  void tearDown() { }

  void doTest(int num);
  void test1() { doTest(1); }
  void test2() { doTest(2); }
  void test3() { doTest(3); }
  void test4() { doTest(4); }
  void test5() { doTest(5); }
  void test6() { doTest(6); }
  void test7() { doTest(7); }
  void test8() { doTest(8); }
  void test9() { doTest(9); }
  void test10() { doTest(10); }
  void test11() { doTest(11); }

  CPPUNIT_TEST_SUITE(TestSegmentsIntersector);
  CPPUNIT_TEST(test1);
  CPPUNIT_TEST(test2);
  CPPUNIT_TEST(test3);
  CPPUNIT_TEST(test4);
  CPPUNIT_TEST(test5);
  CPPUNIT_TEST(test6);
  CPPUNIT_TEST(test7);
  CPPUNIT_TEST(test8);
  CPPUNIT_TEST(test9);
  CPPUNIT_TEST(test10);
  CPPUNIT_TEST(test11);
  CPPUNIT_TEST_SUITE_END();

public:
  static bool verify(const SegmentsIntersector::ResultType &r1,
                     const SegmentsIntersector::ResultType &r2,
                     Real tol);
  static vector<Segment<2>> readInputFile(const string &filename, Real &tol);
  static SegmentsIntersector::ResultType readAnswerFile(const string &filename);
};

#endif // TESTSEGMENTSINTERSECTOR_H
