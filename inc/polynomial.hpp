#pragma once

#include "data_type.hpp"
#include <array>
#include <stdint.h>

namespace math {
namespace polynomial {

template <FloatType T, int32_t N> class Polynomial1D {
public:
  typedef T DateType;

  void eval(const T t) const;

  const T *coeff_address() const { return coeff_.data(); }

private:
  static const int32_t order_ = N;
  static const int32_t coeff_num_ = order_ + 1;

  std::array<T, N + 1> coeff_;
};
} // namespace polynomial
} // namespace math