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
  std::invoke(&T::eval, polynomial, typename T::DateType());
  std::invoke(&T::coeff_address, polynomial);
};

} // namespace polynomial
} // namespace math
