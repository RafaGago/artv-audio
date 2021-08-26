#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// x^(2/3)
struct compand_1a_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return sgn_no_zero (x) * pow (abs (x), (T) (2. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;
    return vec_sgn_no_zero (x) * vec_pow (vec_abs (x), (T) (2. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    //  (3/5)x^(5/3)
    auto ax = abs (x);
    return pow (ax, (T) (2. / 3.)) * ax * (T) (3. / 5.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;
    auto ax = vec_abs (x);
    return vec_pow (ax, (T) (2. / 3.)) * ax * (T) (3. / 5.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    static_assert (sizeof (T) != 0, "sign to be checked!");
    //  (9/40)x^(8/3)
    auto ax = abs (x);
    return pow (ax, (T) (2. / 3.)) * ax * x * (T) (9. / 40.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;
    static_assert (sizeof (T) != 0, "sign to be checked!");
    auto ax = vec_abs (x);
    return vec_pow (ax, (T) (2. / 3.)) * ax * x * (T) (9. / 40.);
  }
  //----------------------------------------------------------------------------
};

// x^(3/2)
struct compand_1b_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return sqrt (abs (x)) * x;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;
    return vec_sqrt (vec_abs (x)) * x;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    // (2/5)x^(5/2)
    auto ax = abs (x);
    return sqrt (ax) * ax * ax * (T) (2. / 5.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;
    auto ax = vec_abs (x);
    return vec_sqrt (ax) * ax * ax * (T) (2. / 5.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    // (4/35)x^(7/2)
    static_assert (sizeof (T) != 0, "sign to be checked!");
    return sqrt (abs (x)) * x * x * x * (T) (4. / 35.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;
    static_assert (sizeof (T) != 0, "sign to be checked!");
    return vec_sqrt (vec_abs (x)) * x * x * x * (T) (4. / 35.);
  }
  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
template <uint order>
using compand_1a_adaa = adaa::waveshaper<compand_1a_functions, order>;

template <uint order>
using compand_1b_adaa = adaa::waveshaper<compand_1b_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
