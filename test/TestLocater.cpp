#include "TestLocater.h"

#include <fstream>
#include <string>
#include <vector>

#include "Core/dirConfig.h"
#include "YinSet/OrientedJordanCurve.h"
#include "YinSet/PointsLocater.h"
#include "YinSet/SegmentedRealizableSpadjor.h"

using std::ifstream;
using std::string;
using std::vector;
using rVec = Vec<Real, 2>;

const Real tol = 1e-8;
// const int Dim = 2;
const int Order = 2;

std::vector<Segment<2>> collapseToSeg(
    const std::vector<OrientedJordanCurve<2, 2>>& bdries,
    std::vector<int>& cIdx, std::vector<int>& pIdx);

void TestLocater::doTest(int num) {
  int numCurve;
  vector<OrientedJordanCurve<2, 2>> jordanCurves;
  vector<rVec> queries;
  vector<int> answer;

  ifstream infile(string(std::string(ROOT_DIR) + "/test/data/testLocater-") +
                  (char)('0' + num) + ".input");
  infile >> numCurve;
  for (int n = 0; n < numCurve; ++n) {
    int nk;
    infile >> nk;
    vector<rVec> knots(nk + 1);
    for (int k = 0; k < nk; ++k) infile >> knots[k][0] >> knots[k][1];
    knots.back() = knots.front();
    jordanCurves.push_back(
        OrientedJordanCurve<2, 2>(fitCurve<2>(knots, Curve<2, 2>::periodic)));
  }
  SegmentedRealizableSpadjor<Order> srs(jordanCurves);

  int nq;
  infile >> nq;
  queries.resize(nq);
  for (int k = 0; k < nq; ++k) infile >> queries[k][0] >> queries[k][1];

  infile.close();
  infile.open(string(std::string(ROOT_DIR) + "/test/data/testLocater-") +
              (char)('0' + num) + ".answer");
  answer.resize(nq);
  for (int k = 0; k < nq; ++k) infile >> answer[k];
  infile.close();

  PointsLocater locater(tol);
  std::vector<int> ci, pi;
  auto result =
      locater(collapseToSeg(jordanCurves, ci, pi), queries, srs.isBounded(tol));
  for (int k = 0; k < nq; ++k) {
    CPPUNIT_ASSERT(answer[k] == result[k]);
  }
}
