#include <cassert>
#include <fstream>
#include <iostream>

#include "YinSet/SegmentedRealizableSpadjor.h"
#include "YinSet/YinSet.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

const int Dim = 2;
const int Order = 2;
const Real tol = 1e-8;

int main(int argc, char* argv[]) {
  YinSet<Dim, Order>*panda, *mickey;
  // ===================================================================
  // Load the Yin sets and calculate the Hasse diagram.
  {
    ifstream input1("data1/Panda.dat", std::ios::binary);
    assert(input1);
    panda = new YinSet<Dim, Order>(input1, tol);
    cout << "Hasse of panda : \n";
    cout << panda->getHasseString() << endl;
    cout << "Betti numbers = " << panda->getBettiNumber(0) << ", "
         << panda->getBettiNumber(1) << endl
         << endl;
  }
  {
    ifstream input2("data/mickey.input.dat", std::ios::binary);
    assert(input2);
    mickey = new YinSet<Dim, Order>(input2, tol);
    //    *mickey = mickey->translate({-1.0, 0.0});
    cout << "Hasse of mickey : \n";
    cout << mickey->getHasseString() << endl;
    cout << "Betti numbers = " << mickey->getBettiNumber(0) << ", "
         << mickey->getBettiNumber(1) << endl
         << endl;
  }

  // ===================================================================
  // Unary Op
  //  {
  //    auto panda_comp = panda.complement();
  //    ofstream output("results/panda_complement.dat", std::ios::binary);
  //    panda_comp.save(output);
  //  }
  //  {
  //    auto mickey_comp = mickey.complement();
  //    ofstream output("results/mickey_complement.dat", std::ios::binary);
  //    mickey_comp.save(output);
  //  }

  // ===================================================================
  // Test intersection
  {
    auto cap = intersect(*panda, *mickey, tol);
    cout << "Hasse of panda-mickey intersection: \n";
    cout << cap.getHasseString() << endl;
    cout << "Betti numbers = " << cap.getBettiNumber(0) << ", "
         << cap.getBettiNumber(1) << endl
         << endl;

    ofstream output("results/panda-mickey-ints.dat", std::ios::binary);
    cap.dump(output);
  }

  // ===================================================================
  // Test union
  {
    SegmentedRealizableSpadjor<Order> sPanda(panda->getBoundaryCycles());
    SegmentedRealizableSpadjor<Order> sMickey(mickey->getBoundaryCycles());
    auto sCup =
        meet(sPanda.complement(), sMickey.complement(), tol).complement();
    YinSet<Dim, Order> cup(sCup, tol);
    cout << "Hasse of panda-mickey union:\n";
    cout << cup.getHasseString() << endl;
    cout << "Betti numbers = " << cup.getBettiNumber(0) << ", "
         << cup.getBettiNumber(1) << endl
         << endl;

    ofstream output("results/panda-mickey-union.dat", std::ios::binary);
    cup.dump(output);
  }

  delete mickey;
  delete panda;

  return 0;
}
