#pragma once

#include "data_type.hpp"
#include "eigen3/unsupported/Eigen/Polynomials"
#include <optional>
#include <stdint.h>
#include <limits>
#include <iostream>

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
    if (std::fabs(coeff_ptr[Deg]) <= std::numeric_limits<F>::min()) {
      return PolynomialRootFinder<F, Deg - 1>::FindAllRealRoot(polynomial, root_num, real_roots);
    } else {
      //todo: using eigen get roots
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
    if (std::fabs(coeff_ptr[1]) <= std::numeric_limits<F>::min()) {
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
    assert(polynomial.order() >= Deg);
    const F* coeff_ptr = polynomial.coeff_address();
    if (std::fabs(coeff_ptr[2]) <= std::numeric_limits<F>::min()) {
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


template <FloatType F> class PolynomialRootFinder<F, 3> {
  public:
  static std::optional<RootSolverStatus>
  FindAllRealRoot(const PolynomialItf auto &polynomial, int32_t *const root_num, 
                  F * const real_roots) {
    assert(polynomial.order() >= Deg);
    const F* coeff_ptr = polynomial.coeff_address();
    if (std::fabs(coeff_ptr[3]) <= std::numeric_limits<F>::min()) {
      return PolynomialRootFinder<F, 2>::FindAllRealRoot(polynomial, root_num, real_roots);
    } else {
      const F a = coeff_ptr[3];
      
      const F b = coeff_ptr[2] / a; // speed up
      const F b_sq = b * b;
      const F b_cub = b_sq * b;

      const F c = coeff_ptr[1] / a;

      const F d = coeff_ptr[0] / a;

      const F bc = b * c; 

      const F alpha = bc / 6 - b_cub / 27 - d / 2;
      const F alpha_sq = alpha * alpha;
      const F beta = c / 3 - b_sq / 9;
      const F beta_cub = beta * beta * beta;
      const F delta = alpha_sq + beta_cub;
      const F base = - b / 3;
      if (delta > std::numeric_limits<F>::min()) {
        // delta is positive: one roots
        *root_num = 1;
        
        const F delta_sqrt = std::sqrt(delta);
        real_roots[0] = base + std::cbrt(alpha + delta_sqrt) + std::cbrt(alpha - delta_sqrt);

      } else if (delta < - std::numeric_limits<F>::min()){
        // delta is negative: three roots
        *root_num = 3;

        const F beta_abs = std::fabs(beta);
        const F beta_abs_sqrt = std::sqrt(beta_abs);
        const F beta_abs_sqrt_cub = beta_abs_sqrt * beta_abs_sqrt * beta_abs_sqrt;
        const F theta = std::acos(alpha / beta_abs_sqrt_cub);

        const F scale = beta_abs_sqrt + beta_abs_sqrt;
        
        static constexpr F k2Pi = M_PI + M_PI;
        real_roots[0] = base + scale * std::cos(theta / 3);
        real_roots[1] = base + scale * std::cos((theta + k2Pi) / 3);
        real_roots[2] = base + scale * std::cos((theta - k2Pi) / 3);
      } else  {
        // delta is zero: three same roots or one + two same roots
        const F base = - b / 3;
        const F alpha_cbrt = std::cbrt(alpha); 
        real_roots[0] = base + alpha_cbrt * 2;
        if (alpha < std::numeric_limits<F>::min()) {
          *root_num = 1;
        } else {
          *root_num = 2;
          real_roots[1] = base - alpha_cbrt;
        }
      }
      return RootSolverStatus::kSuccess;
    }
  }
};



} // namespace polynomial
} // namespace math