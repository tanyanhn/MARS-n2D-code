#include <fstream>
#include <sstream>
#include <iomanip>
#include <cctype>
#include "TestSegmentsIntersector.H"
#include "Core/dirConfig.h"

using std::ifstream;
using std::ostringstream;
using std::istringstream;

void TestSegmentsIntersector::doTest(int num)
{
  ostringstream prefix;
  prefix << std::string(ROOT_DIR) + "/test/data/testLineInts-" << num;
  
  Real tol;
  vector<Segment<2>> segs = readInputFile(prefix.str() + ".input", tol);
  auto answer = readAnswerFile(prefix.str() + ".answer");
  SegmentsIntersector si(tol);
  auto output = si(segs);
  CPPUNIT_ASSERT(verify(answer, output, tol));
}

bool my_isspace(const char &c) {
  return isspace(c);
}

vector<Segment<2>> TestSegmentsIntersector::readInputFile(const string &filename, Real &tol)
{
  std::cout << "filename: \n" << filename;
  ifstream input(filename);
  assert(input);
  input >> tol;
  input.ignore(999, '\n'); // or else the first getline returns nothing
  
  vector<Segment<2>> segments;
  string buf;
  Vec<Real,2> p1,p2;
  for(;;) {
    std::getline(input, buf);
    if(!input || all_of(buf.begin(), buf.end(), my_isspace))
      break;
    istringstream iss(buf);
    iss >> p1[0] >> p1[1] >> p2[0] >> p2[1];
    segments.push_back(Segment<2>(p1,p2));
  }
  return segments;
}

SegmentsIntersector::ResultType TestSegmentsIntersector::readAnswerFile(const string &filename)
{
  ifstream input(filename);
  assert(input);

  SegmentsIntersector::ResultType r;
  string buf;
  while(true) {
    getline(input, buf);
    if(!input || all_of(buf.begin(), buf.end(), my_isspace))
      break;
    istringstream iss(buf);
    Vec<Real,2> p;
    vector<int> inc;
    iss >> p[0] >> p[1];
    while(iss) {
      int i = -1;
      iss >> i;
      if(i == -1) break;
      inc.push_back(i);
    }
    r.insert(make_pair(p, inc));
  }
  return r;
}

bool TestSegmentsIntersector::verify(const SegmentsIntersector::ResultType &r1,
                                     const SegmentsIntersector::ResultType &r2,
                                     Real tol)
{
  if(r1.size() != r2.size())
    return false;
  VecCompare<Real,2> pt_cmp(tol);
  auto i = r1.begin();
  for(auto j=r2.begin(); j!=r2.end(); i++, j++) {
    if(pt_cmp.compare(i->first, j->first) != 0)
      return false;
    vector<int> inc1(i->second);
    vector<int> inc2(j->second);
    std::sort(inc1.begin(), inc1.end());
    std::sort(inc2.begin(), inc2.end());
    if(!std::includes(inc1.begin(), inc1.end(), inc2.begin(), inc2.end()))
      return false;
    if(!std::includes(inc2.begin(), inc2.end(), inc1.begin(), inc1.end()))
      return false;
  }
  return true;
}
