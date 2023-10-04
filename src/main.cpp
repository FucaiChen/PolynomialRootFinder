#include "inc/polynomial.hpp"
#include "inc/polynomial_root_finder.hpp"
#include <iostream>

#include <sys/time.h>

using namespace std;

void test_func() {
  typedef double D;
  static constexpr int32_t N = 3;
  std::array<D, N + 1> coeff{0, -1, 0, 1};
  math::polynomial::Polynomial1D<D, N> polynomial1(coeff.data());
  math::polynomial::PolynomialRootFinder<D, N> solver;
  std::int32_t roots_num;
  std::array<D, N> real_roots;

  struct timeval t_begin, t_end;
  static constexpr int32_t kLoop = 1e9;
  gettimeofday(&t_begin,NULL);
  for(int i = 0; i < kLoop; ++i) {
    solver.FindAllRealRoot(polynomial1, &roots_num, &real_roots[0]);
  }
  gettimeofday(&t_end,NULL);

  const double time_us = (t_end.tv_sec - t_begin.tv_sec) * 1e6 +
  (t_end.tv_usec - t_begin.tv_usec);
  std::cout << "time = " << time_us << " us"
  << " av time = " << time_us / kLoop * 1e3 << " ns" << std::endl;

  std::cout << "poly roots num = " << roots_num << std::endl;
  for (int32_t i = 0 ; i < roots_num; ++i) {
    std::cout << "root[" << i << "] = " << real_roots[i] << std::endl;
  }
}

int main(int argc, char **argv) {

  cout << "hello world!" << endl;

  std::cout << "double epsilon = " << std::numeric_limits<double>::epsilon() << std::endl;
  std::cout << "float epsilon = " << std::numeric_limits<float>::epsilon() << std::endl;

  test_func();
  return 0;
}
