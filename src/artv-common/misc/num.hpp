#pragma once

#include <type_traits>

#include "artv-common/misc/vec_util.hpp"

namespace artv {
//------------------------------------------------------------------------------
// A number wrapper.
//
// Sometimes there are numeric classes where overloading the operators to take
// builtin scalar types directly is not desired to avoid errors.
//
// This class is a wrapper around those classes to signal the intent that it is
// desired to operate with a wrapper but not the type directly.
template <class T>
struct num {
  T value;
};

#if defined(__cpp_deduction_guides)
template <class T>
num (T) -> num<T>;
#else
#error "This span implementation requires __cpp_deduction_guides"
#endif

//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator+ (num<T> a, U b)
{
  return vec_cast<U> (a.value) + b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator+ (U a, num<T> b)
{
  return a + vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator- (num<T> a, U b)
{
  return vec_cast<U> (a.value) - b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator- (U a, num<T> b)
{
  return a - vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator* (num<T> a, U b)
{
  return vec_cast<U> (a.value) * b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator* (U a, num<T> b)
{
  return a * vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator/ (num<T> a, U b)
{
  return vec_cast<U> (a.value) / b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator/ (U a, num<T> b)
{
  return a / vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator== (num<T> a, U b)
{
  return vec_cast<U> (a.value) == b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator== (U a, num<T> b)
{
  return a == vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator!= (num<T> a, U b)
{
  return vec_cast<U> (a.value) != b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator!= (U a, num<T> b)
{
  return a != vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator<= (num<T> a, U b)
{
  return vec_cast<U> (a.value) <= b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator<= (U a, num<T> b)
{
  return a <= vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator<(num<T> a, U b)
{
  return vec_cast<U> (a.value) < b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator<(U a, num<T> b)
{
  return a < vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator>= (num<T> a, U b)
{
  return vec_cast<U> (a.value) >= b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator>= (U a, num<T> b)
{
  return a >= vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator> (num<T> a, U b)
{
  return vec_cast<U> (a.value) > b;
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_vec_or_scalar_v<U>>* = nullptr>
constexpr inline U operator> (U a, num<T> b)
{
  return a > vec_cast<U> (b.value);
}
//------------------------------------------------------------------------------
} // namespace artv
