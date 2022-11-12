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
  using int_type                = std::intmax_t;
  using uint_type               = std::uintmax_t;
  using std_ratio               = typename std::ratio<Num, Den>::type;
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
  // will this need to be explicit to avoid ambiguity?
  constexpr operator float()
  {
    return (float) std_ratio::num / (float) std_ratio::den;
  }
  //----------------------------------------------------------------------------
  constexpr operator double()
  {
    return (double) std_ratio::num / (double) std_ratio::den;
  }
  //----------------------------------------------------------------------------
  constexpr operator long double()
  {
    return (long double) std_ratio::num / (long double) std_ratio::den;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <class T>
struct is_artv_ratio : public std::false_type {};

template <std::intmax_t N, std::intmax_t D>
struct is_artv_ratio<ratio<N, D>> : public std::true_type {};

template <class T>
static constexpr bool is_artv_ratio_v = is_artv_ratio<T>::value;
//------------------------------------------------------------------------------
template <class T>
struct is_std_ratio : public std::false_type {};

template <std::intmax_t N, std::intmax_t D>
struct is_std_ratio<std::ratio<N, D>> : public std::true_type {};

template <class T>
static constexpr bool is_std_ratio_v = is_std_ratio<T>::value;
//------------------------------------------------------------------------------
} // namespace artv
