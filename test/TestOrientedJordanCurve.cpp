#include "TestOrientedJordanCurve.h"
#include <string>
#include "YinSet/OrientedJordanCurve.h"

using std::string;

void TestOrientedJordanCurve::doTest(int num) {
  string c1 = "4 1 0 0 1 -1 0 0 -1";
  OrientedJordanCurve<2, 2> curve1;
  curve1.define(c1);
  string c2 = "4 1 0 0 1 -1 0 0 -1";
  OrientedJordanCurve<2, 4> curve2;
  curve2.define(c2);
  OrientedJordanCurve<2, 4>* circlePtr;
  circlePtr = new Circle<4>();
  string c3 = "2 3 2 1 0.1";
  circlePtr->define(c3);
  OrientedJordanCurve<2, 2>* rectanglePtr;
  rectanglePtr = new Rectangle<2>();
  string c4 = "-3 -5 8 11 1 1 0.1";
  rectanglePtr->define(c4);
};

int main() {
  TestOrientedJordanCurve test1;
  test1.doTest(0);
}