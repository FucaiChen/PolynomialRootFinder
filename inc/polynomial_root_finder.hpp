#pragma once

#include "data_type.hpp"
#include <optional>
#include <stdint.h>

namespace math {
namespace polynomial {

template <FloatType F, int32_t N> class PolynomialRootFinder {
public:
  enum RootStatus {
    kErrorInput = 1,
    kSuccess = 2,
  };

  // std::optional<RootStatus>
  void FindAllRealRoot(const PolynomialItf auto &polynomial) {
    const auto coeff_ptr = polynomial.coeff_address();

    return;
  };
};

} // namespace polynomial
} // namespace math