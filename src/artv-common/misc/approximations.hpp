#pragma once

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv::taylor {
//------------------------------------------------------------------------------
// sin with optimization point being zero
// see
// https://www.kvraudio.com/forum/viewtopic.php?p=7338746&sid=234f68ea04eb937b3e0d2652cd72b1e6#p7338746
struct sin_zero {
  //----------------------------------------------------------------------------
  static constexpr bool has_even = false; // has even orders (0 and 1 excluded)
  static constexpr bool has_odd  = true; // has odd orders (0 and 1 excluded)
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr auto table()
  {
    return make_array<T> (
      (T) 0., // order 0 -> k
      (T) 1., // order 1
      (T) -0.166666666666666666666666666667, // order 3...
      (T) 0.00833333333333333333333333333333, // order 5...
      (T) -0.000198412698412698412698412698413,
      (T) 0.0000027557319223985890652557319224,
      (T) -0.0000000250521083854417187750521083854,
      (T) 1.60590438368216145993923771702e-10,
      (T) -7.64716373181981647590113198579e-13,
      (T) 2.81145725434552076319894558301e-15,
      (T) -8.22063524662432971695598123687e-18,
      (T) 1.95729410633912612308475743735e-20,
      (T) -3.86817017063068403771691193152e-23,
      (T) 6.4469502843844733961948532192e-26,
      (T) -9.18368986379554614842571683647e-29,
      (T) 1.13099628864477169315587645769e-31);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// sinh with the optimization point being zero
struct sinh_zero {
  //----------------------------------------------------------------------------
  static constexpr bool has_even = false; // has even orders (0 and 1 excluded)
  static constexpr bool has_odd  = true; // has odd orders (0 and 1 excluded)
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr auto table()
  {
    return make_array<T> (
      (T) 0., // order 0 -> k
      (T) 1., // order 1
      (T) 0.166666666666666666666666666667, // order 3...
      (T) 0.00833333333333333333333333333333, // order 5...
      (T) 0.000198412698412698412698412698413,
      (T) 0.0000027557319223985890652557319224,
      (T) 0.0000000250521083854417187750521083854,
      (T) 1.60590438368216145993923771702e-10,
      (T) 7.64716373181981647590113198579e-13,
      (T) 2.81145725434552076319894558301e-15,
      (T) 8.22063524662432971695598123687e-18,
      (T) 1.95729410633912612308475743735e-20,
      (T) 3.86817017063068403771691193152e-23,
      (T) 6.4469502843844733961948532192e-26,
      (T) 9.18368986379554614842571683647e-29,
      (T) 1.13099628864477169315587645769e-31);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// convert from order to number of terms (of order 2 and above)
template <class Table>
static constexpr uint order_to_n_terms (uint order)
{
  if (order <= 1) {
    return 0;
  }
  return (order - Table::has_odd)
    / ((Table::has_even && Table::has_odd) ? 1 : 2);
}
//------------------------------------------------------------------------------
// N_terms = number of terms of order 2 or bigger that have a nonzero
// coefficient.
//----------------------------------------------------------------------------
template <
  class Table,
  uint N_terms,
  class V,
  enable_if_vec_of_float_point_t<V>* = nullptr>
static V get (V x)
{
  using T = vec_value_type_t<V>;
  static_assert (Table::has_even || Table::has_odd);
  constexpr auto table = Table::template table<T>();
  static_assert (table.size() >= 2);
  static_assert (N_terms <= (table.size() - 2));

  V ret = table[0] + table[1] * x;
  if constexpr (table.has_even && table.has_odd) {
    V pwr = x; // first term in the loop is order 2
    for (uint i = 0; i < N_terms; ++i) {
      pwr *= x;
      ret += table[i] * pwr;
    }
    return ret;
  }
  else {
    V pwr = table.has_even ? vec_set<V> (1) : x;
    for (uint i = 0; i < N_terms; ++i) {
      pwr *= x * x;
      ret += table[i] * pwr;
    }
    return ret;
  }
}
//------------------------------------------------------------------------------

} // namespace artv::taylor
