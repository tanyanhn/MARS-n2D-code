#include <fstream>
#include "TestYinSetsBooleanOps.H"
#include "YinSet/SegmentedRealizableSpadjor.h"
#include "YinSet/YinSet.h"
using std::ifstream;
using std::ofstream;

void TestYinSetsBooleanOps::doTest(const string &name1, const string &name2, Real tol,
                                   const rVec &ofs, const string &answer)
{
  ifstream input1(name1, std::ios::binary);
  ifstream input2(name2, std::ios::binary);
  ifstream input3(answer, std::ios::binary);
  assert(input1 && input2 && input3);

  const int Dim = 2;
  const int Order = 2;
  SegmentedRealizableSpadjor<Order> srs1(input1);
  SegmentedRealizableSpadjor<Order> srs2(input2);
  if(norm(ofs, 0) != 0.0)
    srs2 = srs2.translate(ofs);
  YinSet<Dim, Order> ys3(input3, tol);
  YinSet<Dim, Order> ys4(meet(srs1, srs2, tol), tol);
  if(get_dbglevel() >= 2) {
    ofstream of("results/resultYinSet.dat", std::ios_base::binary);
    ys4.dump(of);
  }
  CPPUNIT_ASSERT(ys3.equal(ys4, tol));
}
