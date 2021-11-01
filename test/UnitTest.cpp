#include <iostream>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCaller.h>
#include "TestSegmentsIntersector.H"
#include "TestLocater.h"
#include "TestYinSetsBooleanOps.H"

int main(int argc, char *argv[])
{
  CppUnit::TextUi::TestRunner runner1;
  runner1.addTest(TestSegmentsIntersector::suite());
  CppUnit::TextUi::TestRunner runner3;
  runner3.addTest(TestLocater::suite());
  CppUnit::TextUi::TestRunner runner2;
  runner2.addTest(TestYinSetsBooleanOps::suite());

  int result = 0;

  std::cout << "TestSegmentsIntersector" << std::endl;
  result -= !runner1.run();
  std::cout << "TestLocater" << std::endl;
  result -= !runner3.run();
  std::cout << "TestYinSetsBooleanOps" << std::endl;
  result -= !runner2.run();

  return result;
}
