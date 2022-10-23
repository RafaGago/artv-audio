#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
// sqrt waveshaper with first derivative AA. Based on this thread:
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=521377&sid=18fa45082ea11269f4c29c3ccf36becd
struct sqrt_sigmoid_functions {
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;

    return x / vec_sqrt (x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_sqrt (x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return ((x * vec_sqrt (x * x + T (1.))) + vec_asinh (x)) * T (0.5);
  }
};
//------------------------------------------------------------------------------
// This is special cased, so it doesn't have to branch on the problematic
// division:
// https://www.kvraudio.com/forum/viewtopic.php?t=521377&start=30
// message by "martinvicanek"
// (sqrt(1 + x^2) - sqrt(1+ x1^2))/(x - x1) =
//    (x + x1)/(sqrt(1 + x^2) + sqrt(1 + x1^2))

class sqrt_sigmoid_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { x1, x1_sqrt, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V tick (xspan<V const>, xspan<V> st, V x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= n_states);

    V x1v     = st[x1];
    V x1vsqrt = st[x1_sqrt];

    V xsqrt = vec_sqrt (x * x + (T) 1.);

    st[x1]      = x;
    st[x1_sqrt] = xsqrt;

    return (x + x1v) / (xsqrt + x1vsqrt);
  }
};

#if 0
template <uint order>
using sqrt_sigmoid_adaa = std::conditional_t<
  order == 1,
  sqrt_sigmoid_adaa_1,
  adaa::waveshaper<sqrt_sigmoid_functions, order>>;
#else
template <uint order>
using sqrt_sigmoid_adaa = adaa::waveshaper<sqrt_sigmoid_functions, order>;
#endif
//------------------------------------------------------------------------------
} // namespace artv
