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

enum FValueSign {
  kPositive = 1,
  kNegative = 2,
  kZero = 3,
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
      //todo: using Eigen get roots
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
      const F den = a + a;
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

    static constexpr F k2Pi = M_PI + M_PI;
    static constexpr F kPiFrac23 = k2Pi / 3;
    static constexpr F kCosPiFrac23 = std::cos(kPiFrac23);
    static constexpr F kSinPiFrac23 = std::sin(kPiFrac23);

    const F* coeff_ptr = polynomial.coeff_address();
    if (std::fabs(coeff_ptr[3]) <= std::numeric_limits<F>::epsilon()) {
      return PolynomialRootFinder<F, 2>::FindAllRealRoot(polynomial, root_num, real_roots);
    } else {
      const F& a = coeff_ptr[3];
      
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

      if (delta > std::numeric_limits<F>::epsilon()) {
        // delta is positive: one roots
        *root_num = 1;
        
        const F delta_sqrt = std::sqrt(delta);
        real_roots[0] = base + std::cbrt(alpha + delta_sqrt) + std::cbrt(alpha - delta_sqrt);

      } else if (delta < - std::numeric_limits<F>::epsilon()){
        // delta is negative: three roots
        *root_num = 3;

        const F beta_abs = std::fabs(beta);
        const F beta_abs_sqrt = std::sqrt(beta_abs);
        const F beta_abs_sqrt_cub = beta_abs_sqrt * beta_abs_sqrt * beta_abs_sqrt;
        const F theta = std::acos(alpha / beta_abs_sqrt_cub);

        const F scale = beta_abs_sqrt + beta_abs_sqrt;
        
        const F gamma = theta / 3;
        const F s_cos_gamma = scale * std::cos(gamma);
        const F s_sin_gamma = scale * std::sin(gamma);

        const F cc = s_cos_gamma * kCosPiFrac23;
        const F ss = s_sin_gamma * kSinPiFrac23;
        real_roots[0] = base + s_cos_gamma;//base + scale * std::cos(theta / 3);
        real_roots[1] = base + (cc + ss);// base + scale * std::cos((theta + k2Pi) / 3);
        real_roots[2] = base + (cc - ss);// base + scale * std::cos((theta - k2Pi) / 3);
      } else  {
        // delta is zero: three same roots or one + two same roots
        const F alpha_cbrt = std::cbrt(alpha); 
        real_roots[0] = base + alpha_cbrt * 2;
        if (alpha < std::numeric_limits<F>::epsilon()) {
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



template <FloatType F> class PolynomialRootFinder<F, 4> { 
  public:
  static std::optional<RootSolverStatus>
  FindAllRealRoot(const PolynomialItf auto &polynomial, int32_t *const root_num, 
                  F * const real_roots) {
    assert(polynomial.order() >= Deg);

    const F* coeff_ptr = polynomial.coeff_address();
    if (std::fabs(coeff_ptr[4]) <= std::numeric_limits<F>::epsilon()) {
      return PolynomialRootFinder<F, 3>::FindAllRealRoot(polynomial, root_num, real_roots);
    } else {
      const F& a = coeff_ptr[4];
      const F b = coeff_ptr[3] / a;
      const F c = coeff_ptr[2] / a;
      const F d = coeff_ptr[1] / a; 
      const F e = coeff_ptr[0] / a; 

      const F b_sq = b * b;
      const F b_cub = b_sq * b;
      const F b_qua = b_sq * b_sq;

      const F bc = b * c;
      const F b_sq_c = b_sq * c;
      
      const F bd = b * d;

      static constexpr F n_frac_3_8 = - static_cast<F>(3) / static_cast<F>(8);
      static constexpr F n_3_256 = - static_cast<F>(3) / static_cast<F>(256);
      
      const F alpha = n_frac_3_8 * b_sq + c;
      const F alpha_sq = alpha * alpha;
      const F alpha_cub = alpha_sq * alpha;

      const F beta = b_cub / 8 - bc / 2 + d;
      const F gamma = n_3_256 * b_qua + b_sq_c / 16 - bd / 4 + e;

      const F base = - b / 4;

      if (std::fabs(beta) < std::numeric_limits<F>::epsilon()) { 
        // beta == 0
        const F alpha_sq_M_4_gamma = alpha_sq - 4 * gamma;
        
        if (alpha_sq_M_4_gamma < - std::numeric_limits<F>::epsilon()) {
          // exist none root
          *root_num = 0;
          return RootSolverStatus::kErrorInput;
        } else if (alpha_sq_M_4_gamma > std::numeric_limits<F>::epsilon()){
          *root_num = 2;
          const F offset = std::sqrt(std::fabs(alpha - std::sqrt(alpha_sq_M_4_gamma))/ 2);
          real_roots[0] = base + offset;
          real_roots[1] = base - offset;
          return RootSolverStatus::kSuccess;
        } else {
          *root_num = 1;
          real_roots[0] = base;
          return RootSolverStatus::kSuccess;
        }

      } else {
        // beta != 0
        const F P = - alpha_sq / 12 - gamma;
        const F P_cub = P * P * P;
        const F Q = - alpha_cub / 108 + alpha * gamma / 3 - beta * beta / 8;
        const F Q_sq = Q * Q;
        F temp = Q_sq / 4 + P_cub / 27;
        if (temp < - std::numeric_limits<F>::epsilon()) {
          *root_num = 0;
          return RootSolverStatus::kErrorInput;
        }
        temp = std::max(F(), temp);
        const F R = Q / 2 + std::sqrt(temp);
        
        if (R < -std::numeric_limits<F>::epsilon()) {
          *root_num = 0;
          return RootSolverStatus::kErrorInput;
        }
        const F U = std::cbrt(R);

        static constexpr F n_frac_5_6 = - static_cast<F>(5) / static_cast<F>(6);
        F y = n_frac_5_6 * alpha;
        if (std::fabs(U) < std::numeric_limits<F>::epsilon()) {
          y -= std::cbrt(Q);
        } else {
          y += U - P / U / 3;
        }

        const F alpha_2y = alpha + y + y;
        if (alpha_2y < std::numeric_limits<F>::epsilon()) {
          *root_num = 0;
          return RootSolverStatus::kErrorInput;
        }
        const F alpha_2y_sqrt = std::sqrt(alpha_2y);
        
        const F alpha3_2y = alpha_2y + alpha + alpha;
        const F beta2_div_alpha_2y_sqrt = (beta + beta) / alpha_2y_sqrt;
        
        const F temp_p = - (alpha3_2y + beta2_div_alpha_2y_sqrt);
        const F temp_n = - (alpha3_2y - beta2_div_alpha_2y_sqrt);

        *root_num = 0;
        if (temp_p > std::numeric_limits<F>::epsilon()) {
          const F temp_p_sqrt = std::sqrt(temp_p);
          real_roots[*root_num] = base + (alpha_2y_sqrt + temp_p_sqrt) / 2;
          ++(*root_num);
          real_roots[*root_num] = base + (alpha_2y_sqrt - temp_p_sqrt) / 2;
          ++(*root_num);
        } else if (temp_p < - std::numeric_limits<F>::epsilon()) {
          // none root for this case
        } else {
          real_roots[*root_num] = base + alpha_2y_sqrt / 2;
          ++(*root_num);
        }

        if (temp_n > std::numeric_limits<F>::epsilon()) {
          const F temp_n_sqrt = std::sqrt(temp_n);
          real_roots[*root_num] = base + (- alpha_2y_sqrt + temp_n_sqrt) / 2;
          ++(*root_num);
          real_roots[*root_num] = base + (- alpha_2y_sqrt - temp_n_sqrt) / 2;
          ++(*root_num);
        } else if (temp_p < - std::numeric_limits<F>::epsilon()) {
          // none root for this case
        } else {
          if (std::fabs(temp_p) > std::numeric_limits<F>::epsilon()) {
            real_roots[*root_num] = base + alpha_2y_sqrt / 2;
            ++(*root_num);
          }
        }
        return RootSolverStatus::kSuccess;
      }
    }

  }
};


} // namespace polynomial
} // namespace math