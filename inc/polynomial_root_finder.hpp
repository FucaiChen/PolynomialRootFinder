#pragma once

#include "data_type.hpp"
#include "eigen3/unsupported/Eigen/Polynomials"
#include <optional>
#include <stdint.h>
#include <limits>

namespace math {
namespace polynomial {

enum RootSolverStatus {
  kErrorInput = 1,
  kSuccess = 2,
};

template <FloatType F, int32_t Deg> class PolynomialRootFinder {
public:

  /*/**
   * @brief find all the real root of a given polynomial
   * 
   * @param polynomial 
   * @return std::optional<RootSolverStatus> 
   */
  static std::optional<RootSolverStatus>
  FindAllRealRoot(const PolynomialItf auto &polynomial, int32_t *const root_num, F * const real_roots) {
    assert(polynomial.order() >= Deg);
    const auto coeff_ptr = polynomial.coeff_address();
    if (coeff_ptr[Deg] == F()) {
      return PolynomialRootFinder<F, Deg - 1>::FindAllRealRoot(polynomial, root_num, real_roots);
    } else {
      return RootSolverStatus::kSuccess;
    }
  }

private:
  static Eigen::PolynomialSolver<F, Deg> eigen_root_solver_;
};




template <FloatType F> class PolynomialRootFinder<F, 1> {
  public:
  static std::optional<RootSolverStatus>
  FindAllRealRoot(const PolynomialItf auto &polynomial, int32_t *const root_num, F * const real_roots) {
    assert(polynomial.order() >= Deg);
    const auto coeff_ptr = polynomial.coeff_address();
    if (coeff_ptr[1] == F()) {
      *root_num = 0;
      real_roots[0] = std::numeric_limits<F>::infinity();
      return RootSolverStatus::kErrorInput;
    } else {
      *root_num = 1;
      real_roots[0] = - coeff_ptr[0] / coeff_ptr[1];
      return RootSolverStatus::kSuccess;
    }
  }
};

template <FloatType F> class PolynomialRootFinder<F, 2> {
  public:
  static std::optional<RootSolverStatus>
  FindAllRealRoot(const PolynomialItf auto &polynomial, int32_t *const root_num, F * const real_roots) {
    assert(polynomial.order() == Deg);
    const F* coeff_ptr = polynomial.coeff_address();
    if (coeff_ptr[2] == F()) {
      return PolynomialRootFinder<F, 1>::FindAllRealRoot(polynomial, root_num, real_roots);
    } else {
      const F& a = coeff_ptr[2];
      const F& b = coeff_ptr[1];
      const F& c = coeff_ptr[0];
      const F den = 2 * a;
      const F delta_sq = b * b - 4 * a * c;
      if (delta_sq < F()) {
        *root_num = 0;
      } else if(delta_sq == F()) {
        *root_num = 1;
        real_roots[0] = - b / den;
      } else {
        *root_num = 2;
        const F delta = std::sqrt(delta_sq);
        const F base = - b / den;
        const F offset = delta / den;
        real_roots[0] = base - offset; // {-b - sqrt(b^2-4ac)} / (2a)
        real_roots[1] = base + offset; // {-b + sqrt(b^2-4ac)} / (2a)
      }
      return RootSolverStatus::kSuccess;
    }
  }
};


} // namespace polynomial
} // namespace math