#pragma once

#include <algorithm>
#include <limits>
#include <ratio>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <std::intmax_t Num, std::intmax_t Den = 1>
class ratio {
public:
  using int_type  = std::intmax_t;
  using uint_type = std::uintmax_t;
  using std_ratio = typename std::ratio<Num, Den>::type;

  static constexpr int_type num = std_ratio::num;
  static constexpr int_type den = std_ratio::den;
  //----------------------------------------------------------------------------
  // A RAW fixed point representation.
  struct fixpt {
    int_type value; // raw fixed point representation
    uint     n_sign; // sign required (if the value is negative)
    uint     n_int; // number of integer bits on "value"
    uint     n_frac; // number of fractional bits on "value"
  };
  //----------------------------------------------------------------------------
  static constexpr fixpt get_fixpt()
  {
    static_assert (
      num != std::numeric_limits<int_type>::min(),
      "Invalid value. Implementation negates the numerator");

    constexpr bool neg = (num < 0);
    // -num for will work from signed to
    constexpr uint_type unum = neg ? -num : num;
    constexpr uint_type uden = den;

    // Integer math for the integer part followed by long division for the
    // decimal.
    uint_type quo = unum / uden;
    uint_type rem = unum % uden;

    uint n_int  = last_bit_set (quo);
    uint n_frac = 0;

    constexpr uint max_iter = (sizeof (int_type) * 8) - 1 /*sign "bit"*/;
    while ((rem != 0) && ((n_int + n_frac) < max_iter)) {
      quo <<= 1;
      rem <<= 1;
      if (rem >= uden) {
        rem -= uden;
        quo += 1;
      }
      ++n_frac;
    }
    return {neg ? -((int_type) quo) : (int_type) quo, neg, n_int, n_frac};
  }
  //----------------------------------------------------------------------------
};

template <typename T, typename = void>
struct is_ratio : std::false_type {};

template <typename T>
struct is_ratio<T, decltype ((void) T::num, (void) T::den, void())>
  : std::true_type {};

template <class T>
static constexpr bool is_ratio_v = is_ratio<T>::value;

namespace detail {
template <class T, class R>
using enable_if_floatpt_and_ratio
  = std::enable_if_t<std::is_floating_point_v<T> && is_ratio_v<R>>;

template <class T, class Ratio>
constexpr inline T ratio_to_float (Ratio)
{
  return (T) Ratio::num / (T) Ratio::den;
}

} // namespace detail
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator+ (T lhs, Ratio rhs) noexcept
{
  return lhs + detail::ratio_to_float<T> (rhs);
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator+ (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) + rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator- (T lhs, Ratio rhs) noexcept
{
  return lhs - detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator- (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) - rhs;
}
//-----------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator* (T lhs, Ratio rhs) noexcept
{
  return lhs * detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator* (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) * rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator/ (T lhs, Ratio rhs) noexcept
{
  return lhs / detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr auto operator/ (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) / rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator== (T lhs, Ratio rhs) noexcept
{
  return lhs == detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator== (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) == rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator!= (T lhs, Ratio rhs) noexcept
{
  return lhs != detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator!= (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) != rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator> (T lhs, Ratio rhs) noexcept
{
  return lhs > detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator> (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) > rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator>= (T lhs, Ratio rhs) noexcept
{
  return lhs >= detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator>= (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) >= rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator<(T lhs, Ratio rhs) noexcept
{
  return lhs < detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator<(Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) < rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator<= (T lhs, Ratio rhs) noexcept
{
  return lhs <= detail::ratio_to_float<T> (rhs);
  ;
}

template <
  class T,
  class Ratio,
  detail::enable_if_floatpt_and_ratio<T, Ratio>* = nullptr>
constexpr bool operator<= (Ratio lhs, T rhs) noexcept
{
  return detail::ratio_to_float<T> (lhs) <= rhs;
}
//------------------------------------------------------------------------------
} // namespace artv
