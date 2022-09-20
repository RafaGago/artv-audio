#pragma once

#include <array>
#include <type_traits>

#include <boost/hana.hpp>

namespace artv {

namespace hana = boost::hana;

//------------------------------------------------------------------------------
template <
  class S,
  class T,
  T... Vs,
  class Offset = std::integral_constant<T, 1>>
constexpr auto array_of_str_with_numeric_suffix (
  S boost_hana_string,
  std::integer_sequence<T, Vs...>,
  Offset offset = Offset {})
{
  static_assert (sizeof...(Vs) < 100, "");
  if constexpr (sizeof...(Vs) == 1) {
    // no suffix if the size is 1
    return std::array<char const*, 1> {boost_hana_string.c_str()};
  }
  else {
    return std::array<char const*, sizeof...(Vs)> {
      (boost_hana_string
       + boost::hana::string_c<
         '_',
         '0' + ((Vs + offset.value) / 10),
         '0' + ((Vs + offset.value) % 10)>)
        .c_str()...};
  }
}
//------------------------------------------------------------------------------
} // namespace artv
