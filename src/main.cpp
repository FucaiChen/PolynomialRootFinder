#include "inc/polynomial.hpp"
#include "inc/polynomial_root_finder.hpp"
#include <iostream>
using namespace std;

void test_func() {
  static constexpr int32_t N = 4;
  std::array<double, N> coeff{0.0, 1.0, 2, 1.0};
  math::polynomial::Polynomial1D<double, N> polynomial1(coeff.data());
  math::polynomial::PolynomialRootFinder<double, N> solver;
  std::int32_t roots_num;
  std::array<double, N> real_roots;
  solver.FindAllRealRoot(polynomial1, &roots_num, &real_roots[0]);

  std::cout << "poly roots num = " << roots_num << std::endl;
  for (int32_t i = 0 ; i < roots_num; ++i) {
    std::cout << "root[" << i << "] = " << real_roots[i] << std::endl;
  }
}

int main(int argc, char **argv) {

  cout << "hello world!" << endl;
  test_func();
  return 0;
}
