#pragma once

#include "data_type.hpp"
#include <array>
#include <stdint.h>

namespace math {
namespace polynomial {

template <FloatType T, int32_t N> class Polynomial1D {
public:
  typedef T DateType;

  explicit Polynomial1D(const T* coeff) {
    for (int32_t i = 0; i <= N; ++i) {
      coeff_[i] = coeff[i];
    }
  }

  void trans(const T delta) {
    coeff_[0] += delta;
  }

  void eval(const T t) const;

  int32_t size() const { return coeff_num_; }

  int32_t order() const {return order_; }

  const T *coeff_address() const { return coeff_.data(); }

private:
  static const int32_t order_ = N;
  static const int32_t coeff_num_ = order_ + 1;

  std::array<T, N + 1> coeff_;
};
} // namespace polynomial
} // namespace math