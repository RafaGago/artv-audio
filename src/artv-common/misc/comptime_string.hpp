#pragma once

// helpers to e.g. parse user-defined literals. This could probably be done with
// boost::hana too.

#include <type_traits>

#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv::comptime {

//------------------------------------------------------------------------------
template <char C>
using chr = std::integral_constant<char, C>;

template <char... Chars>
using str = mp_list<chr<Chars>...>;
//------------------------------------------------------------------------------
constexpr uint exp10 (uint n)
{
  uint exp10 = 1;
  for (uint i = 0; i < n; ++i) {
    exp10 *= 10;
  }
  return exp10;
}
//------------------------------------------------------------------------------
template <class Str>
constexpr std::intmax_t atoi()
{
  constexpr auto len    = mp11::mp_size<Str>::value;
  std::intmax_t  sum    = 0;
  std::intmax_t  factor = 1;

  static_assert (len > 0, "Empty strings not allowed");
  constexpr bool negative = mp11::mp_first<Str>::value == '-';
  using str = std::conditional_t<negative, mp11::mp_pop_front<Str>, Str>;

  factor = comptime::exp10 (len - 1 - negative);
  mp11::mp_for_each<Str> ([&sum, &factor] (auto v) {
    constexpr char vchar = v.value;
    static_assert (vchar >= '0' && vchar <= '9', "Only digits allowed");
    sum += (vchar - '0') * factor;
    factor /= 10;
  });
  return negative ? -sum : sum;
}
//------------------------------------------------------------------------------
} // namespace artv::comptime
