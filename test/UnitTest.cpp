#include <cppunit/TestCaller.h>
#include <cppunit/ui/text/TestRunner.h>
#include <iostream>
#include "TestLocater.h"
#include "TestOrientedJordanCurve.h"
#include "TestSegmentsIntersector.H"
#include "TestYinSetsBooleanOps.H"
int main(int argc, char* argv[]) {
  CppUnit::TextUi::TestRunner runner1;
  runner1.addTest(TestSegmentsIntersector::suite());
  CppUnit::TextUi::TestRunner runner3;
  runner3.addTest(TestLocater::suite());
  CppUnit::TextUi::TestRunner runner2;
  runner2.addTest(TestYinSetsBooleanOps::suite());
  CppUnit::TextUi::TestRunner runner4;
  runner4.addTest(TestOrientedJordanCurve::suite());

  int result = 0;

  std::cout << "TestSegmentsIntersector" << std::endl;
  result -= !runner1.run();
  std::cout << "TestLocater" << std::endl;
  result -= !runner3.run();
  std::cout << "TestYinSetsBooleanOps" << std::endl;
  result -= !runner2.run();
  std::cout << "TestOrientedJordanCurve" << std::endl;
  result -= !runner4.run();
  return result;
}
