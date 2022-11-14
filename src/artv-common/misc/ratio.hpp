#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <ratio>
#include <type_traits>

#include "artv-common/misc/comptime_string.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "boost/mp11/algorithm.hpp"
#include "boost/mp11/list.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <std::intmax_t Num, std::intmax_t Den = 1>
class ratio {
public:
  using int_type  = std::intmax_t;
  using uint_type = std::uintmax_t;
  using std_ratio = typename std::ratio<Num, Den>::type;

  template <std::intmax_t Num_v, std::intmax_t Den_v>
  using rebind = ratio<Num_v, Den_v>;

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

// ratios can't be deduced because the extra parameters might be non-type
// template parameters, hence the need for a member rebind.
template <typename T, typename = void>
struct ratio_can_rebind : std::false_type {};

template <typename T>
struct ratio_can_rebind<
  T,
  decltype ((void) typename T::template rebind<1, 2> {}, void())>
  : std::true_type {};

template <class T>
static constexpr bool ratio_can_rebind_v = ratio_can_rebind<T>::value;

namespace detail {
template <class T, class R>
using enable_if_floatpt_and_ratio
  = std::enable_if_t<std::is_floating_point_v<T> && is_ratio_v<R>>;

template <class R1, class R2>
using enable_if_2x_ratios = std::enable_if_t<is_ratio_v<R1> && is_ratio_v<R2>>;

template <class T, class Ratio>
constexpr inline T ratio_to_float (Ratio)
{
  return (T) Ratio::num / (T) Ratio::den;
}

} // namespace detail
// -----------------------------------------------------------------------------
// Ratios with floating point types
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
// Ratios
//------------------------------------------------------------------------------
template <
  class Ratio1,
  class Ratio2,
  detail::enable_if_2x_ratios<Ratio1, Ratio2>* = nullptr>
constexpr auto operator+ (Ratio1 lhs, Ratio2 rhs) noexcept
{
  if constexpr (ratio_can_rebind_v<Ratio1>) {
    constexpr auto den = Ratio1::den * Ratio2::den;
    constexpr auto num = Ratio1::num * Ratio2::den + Ratio2::num * Ratio1::den;
    return typename Ratio1::template rebind<num, den> {};
  }
  else {
    static_assert (
      std::is_same_v<Ratio1, void>,
      "Ratio operator only available for ratios that can rebind");
  }
}
//------------------------------------------------------------------------------
template <
  class Ratio1,
  class Ratio2,
  detail::enable_if_2x_ratios<Ratio1, Ratio2>* = nullptr>
constexpr auto operator- (Ratio1 lhs, Ratio2 rhs) noexcept
{
  if constexpr (ratio_can_rebind_v<Ratio1>) {
    constexpr auto den = Ratio1::den * Ratio2::den;
    constexpr auto num = Ratio1::num * Ratio2::den - Ratio2::num * Ratio1::den;
    return typename Ratio1::template rebind<num, den> {};
  }
  else {
    static_assert (
      std::is_same_v<Ratio1, void>,
      "Ratio operator only available for ratios that can rebind");
  }
}
//-----------------------------------------------------------------------------
template <
  class Ratio1,
  class Ratio2,
  detail::enable_if_2x_ratios<Ratio1, Ratio2>* = nullptr>
constexpr auto operator* (Ratio1 lhs, Ratio2 rhs) noexcept
{
  if constexpr (ratio_can_rebind_v<Ratio1>) {
    constexpr auto num = Ratio1::num * Ratio2::num;
    constexpr auto den = Ratio1::den * Ratio2::den;
    return typename Ratio1::template rebind<num, den> {};
  }
  else {
    static_assert (
      std::is_same_v<Ratio1, void>,
      "Ratio operator only available for ratios that can rebind");
  }
}
//------------------------------------------------------------------------------
template <
  class Ratio1,
  class Ratio2,
  detail::enable_if_2x_ratios<Ratio1, Ratio2>* = nullptr>
constexpr auto operator/ (Ratio1 lhs, Ratio2 rhs) noexcept
{
  if constexpr (ratio_can_rebind_v<Ratio1>) {
    constexpr auto num = Ratio1::num * Ratio2::den;
    constexpr auto den = Ratio1::den * Ratio2::num;
    return typename Ratio1::template rebind<num, den> {};
  }
  else {
    static_assert (
      std::is_same_v<Ratio1, void>,
      "Ratio operator only available for ratios that can rebind");
  }
}
//------------------------------------------------------------------------------
template <class Ratio, std::enable_if_t<is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator- (Ratio rhs) noexcept
{
  constexpr auto num = -Ratio::num;
  constexpr auto den = Ratio::den;
  return typename Ratio::template rebind<num, den> {};
}
//------------------------------------------------------------------------------
// Literal for creating integer ratios that can then be combined with the
// operators. E.g. for creating a 1/2 ratio instead of doing:
//
// > ratio<1, 2>()
//
// It is possible to do:
//
// 0.5_r
//
// Which what basically does is to count the decimal places and create a ratio
// of 5/10 -> ratio<5,10>. Don't forget that even if the template parameters are
// 5 and 10, "ratio<5,10>::num" is "1" and "ratio<5,10>::den" is "2"; it runs
// the same compile-time GCD present on "std::ratio".
//
// Integers are also supported and the algorithm for creating float ratios is
// relatively naïve, so for e.g. creating a ratio of 1/3 you might want to do
// instead:
//
// (1_r / 3_r)
//
// Which will result in "ratio<1, 3>" instead of e.g.
// "ratio<333333333, 1000000000>"
//
// Implementing something more complicated than this is hairy, as e.g.
// approximating 0.33333333 to 1/3 would require an epsilon/error tolerance.
// That is not easily achieved with user defined literals.

template <char... Chars>
constexpr auto operator"" _r()
{
  constexpr uint len     = sizeof...(Chars);
  using str              = comptime::str<Chars...>;
  using dot              = comptime::chr<'.'>;
  constexpr auto dot_pos = mp11::mp_find<str, dot>::value;
  if constexpr (dot_pos == len) {
    // a plain integer
    return ratio<comptime::atoi<str>(), 1> {};
  }
  else {
    using num_str      = mp11::mp_erase_c<str, dot_pos, dot_pos + 1>;
    constexpr auto num = comptime::atoi<num_str>();
    constexpr uint den = comptime::exp10 (len - 1 - dot_pos);
    return ratio<num, den> {};
  }
}
//------------------------------------------------------------------------------
} // namespace artv
