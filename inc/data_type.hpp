#pragma once

#include <concepts>
#include <functional>
#include <type_traits>
namespace math {
namespace polynomial {

template <typename T>
concept FloatType = std::is_floating_point_v<T>;

template <typename T>
concept PolynomialItf = requires(T polynomial) {
  polynomial.size();
  polynomial.order();
  polynomial.eval(typename T::DateType());
  polynomial.coeff_address();
};

} // namespace polynomial
} // namespace math
