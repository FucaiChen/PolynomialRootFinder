#include "inc/polynomial.hpp"
#include "inc/polynomial_root_finder.hpp"
#include <iostream>

#include <sys/time.h>
#include <random>
#include <algorithm>

#include <eigen3/Eigen/Dense>

using namespace std;


template <typename T, int32_t Deg>
void GeneratePolyCoeff(const std::array<T, Deg>& roots, 
  std::array<T, Deg+1> &coeff) {
  static_assert(Deg > 0);

  std::fill(coeff.begin(), coeff.end() - 1, T());
  coeff.back() = static_cast<T>(1.0);

  int32_t index_begin = Deg - 1;
  
  for (const auto& root: roots) {
    for (int32_t i = index_begin; i < Deg ; ++i) {
      coeff[i] -= coeff[i+1] * root;
    }
    --index_begin;
  }
};


/**
 * @brief test for specified degree polynomial root finder
 * 
 * @param deg 
 */
template <typename T, int32_t Deg>
void TestPolynomialRootFinder() {
  constexpr int32_t N = Deg + 1;


  std::array<T, Deg> roots_bench;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(1, 2);

  std::generate(roots_bench.begin(), roots_bench.end(), 
  [&dis, &gen]()->T{return dis(gen);});

  std::cout << "benchmark roots = ";
  std::for_each(roots_bench.cbegin(), roots_bench.cend(), 
    [](const T& data) {std::cout << data << ", ";});
  std::cout << std::endl;



  std::array<T, N> coeff;
  GeneratePolyCoeff<T, Deg>(roots_bench, coeff);


  math::polynomial::Polynomial1D<T, Deg> polynomial1(coeff.data());
  math::polynomial::PolynomialRootFinder<T, Deg> solver;
  std::int32_t roots_num;
  std::array<T, Deg> real_roots;

  struct timeval t_begin, t_end;
  static constexpr int32_t kLoop = 1;
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
};

void TestCubicPolynomial() {
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

template <typename T>
void TestGeneratePolyCoeff() {
  std::array<T, 3> roots;
  std::array<T, 4> coeff_bench;
  std::array<T, 4> coeff_ret;

  roots[0] = 0;
  roots[1] = 1;
  roots[2] = 2;
  coeff_bench[0] = 0;
  coeff_bench[1] = 2;
  coeff_bench[2] = -3;
  coeff_bench[3] = 1;

  GeneratePolyCoeff<T, 3>(roots, coeff_ret);

  std::cout << "coeff bench\t = ";
  std::for_each(coeff_bench.cbegin(), coeff_bench.cend(), 
  [](const T& data)->void {std::cout << data << ",\t";});
  std::cout << std::endl;

  std::cout << "coeff ret\t = ";
  std::for_each(coeff_ret.cbegin(), coeff_ret.cend(), 
  [](const T& data)->void {std::cout << data << ",\t";});
  std::cout << std::endl;

}

int main(int argc, char **argv) {

  // std::cout << "double epsilon = " << std::numeric_limits<double>::epsilon() << std::endl;
  // std::cout << "float epsilon = " << std::numeric_limits<float>::epsilon() << std::endl;

  TestPolynomialRootFinder<double, 3>();

  // TestCubicPolynomial();

  // TestGeneratePolyCoeff<double>();
  return 0;
}
