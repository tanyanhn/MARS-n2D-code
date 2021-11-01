#ifndef BAYS_TESTLOCATER_H
#define BAYS_TESTLOCATER_H

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class TestLocater : public CppUnit::TestFixture
{
public:
  void setUp() { }
  void tearDown() { }

  void doTest(int num);
  void test1() { doTest(1); }
  void test2() { doTest(2); }
  void test3() { doTest(3); }
  void test4() { doTest(4); }

  CPPUNIT_TEST_SUITE(TestLocater);
    CPPUNIT_TEST(test1);
    CPPUNIT_TEST(test2);
    CPPUNIT_TEST(test3);
    CPPUNIT_TEST(test4);
  CPPUNIT_TEST_SUITE_END();
};

#endif //BAYS_TESTLOCATER_H
