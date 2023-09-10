#include "inc/polynomial.hpp"
#include "inc/polynomial_root_finder.hpp"
#include <iostream>
using namespace std;

void test_func() {

  math::polynomial::Polynomial1D<double, 3> polynomial1;
  math::polynomial::PolynomialRootFinder<double, 3> solver;
  solver.FindAllRealRoot(polynomial1);
}

int main(int argc, char **argv) {

  cout << "hello world!" << endl;
  return 0;
}
