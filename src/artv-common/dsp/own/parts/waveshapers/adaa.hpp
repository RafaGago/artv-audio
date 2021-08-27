#pragma once

// Implementation heavily "inspired" on:
// https://github.com/jatinchowdhury18/ADAA

#include <cmath>

#include "artv-common/dsp/own/parts/filters/moving_average.hpp"
#include "artv-common/dsp/own/parts/misc/interpolators.hpp"

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace adaa {

//------------------------------------------------------------------------------
static constexpr double epsilon (double)
{
  return 0.0000000001;
}
static constexpr float epsilon (float)
{
  return 0.00001;
}
//------------------------------------------------------------------------------
template <class functions, uint order>
class waveshaper;
//------------------------------------------------------------------------------
template <class functions>
class waveshaper<functions, 0> {
public:
  enum coeffs { n_coeffs };
  enum state { n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (crange<vec_value_type_t<V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V                           x)
  {
    return functions::fn (x);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <class functions>
class waveshaper<functions, 1> {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_int, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (crange<vec_value_type_t<V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V                           x)
  {
    // as this has to calculate both branches, it might not be worth bothering.
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (traits.size * n_states));

    T* x1v_ptr     = &st[x1 * traits.size];
    T* x1v_int_ptr = &st[x1_int * traits.size];

    V x1v     = vec_load<V> (x1v_ptr);
    V x1v_int = vec_load<V> (x1v_int_ptr);

    V x_int = functions::int_fn (x);

    vec_store (x1v_ptr, x);
    vec_store (x1v_int_ptr, x_int);

    V diff = x - x1v;
    V big, small;

    auto sel       = vec_abs (diff) >= epsilon (T {});
    bool all_big   = vec_is_all_ones (sel);
    bool all_small = vec_is_all_zeros (sel);

    if (!all_small) {
      big = (x_int - x1v_int) / diff;
    }
    if (!all_big) {
      small = functions::fn ((x + x1v) * (T) 0.5);
    }
    return sel ? big : small;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <class functions>
class waveshaper<functions, 2> {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x2, x2_der, x1_int2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (crange<vec_value_type_t<V>>)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  // kept as an impl reference if going for it.
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    using fns = functions;

    T x2v     = st[x2];
    T x2v_der = st[x2_der];
    T x1v_der = get_derivative_x1 (st, x);
    T diff    = x - x2v;
    T ret;

    if (abs (diff) >= epsilon (T {})) {
      ret = ((T) 2. * (x1v_der - x2v_der)) / diff;
    }
    else {
      // fallback
      T a  = (x + x2v) * (T) 0.5;
      diff = a - x;

      if (abs (diff) >= epsilon (T {})) {
        ret = ((T) 2. / diff)
          * (fns::int_fn (a) + (fns::int2_fn (x) - fns::int2_fn (a)) / diff);
      }
      else {
        ret = fns::fn ((a + x) * (T) 0.5);
      }
    }

    st[x2]     = st[x1];
    st[x1]     = x;
    st[x2_der] = x1v_der;

    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V                           x)
  {
    // Real SIMD: TODO. Might not be worth because of the high number of
    // branches.
    static_assert (!std::is_same_v<V, V>, "To be implemented");
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  static T get_derivative_x1 (crange<T> st, T x)
  {
    T x1v      = st[x1];
    T x1v_int2 = st[x1_int2];

    T x_int2 = functions::int2_fn (x);
    T diff   = x - x1v;

    st[x1_int2] = x_int2;

    if (abs (diff) >= epsilon (T {})) {
      return (x_int2 - x1v_int2) / diff;
    }
    else {
      return functions::int_fn ((x + x1v) * (T) 0.5);
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order, template <uint> class Impl>
class fix_eq_and_delay;

template <template <uint> class Impl>
class fix_eq_and_delay<1, Impl> {
public:
  enum coeffs {
    allpass_coeffs_idx,
    boxcar_coeffs_idx = allpass_coeffs_idx + allpass_interpolator::n_coeffs,
    impl_coeffs_idx   = boxcar_coeffs_idx + moving_average<2>::n_coeffs,
    n_coeffs          = impl_coeffs_idx + Impl<1>::n_coeffs
  };
  enum state {
    allpass_states_idx,
    boxcar_states_idx = allpass_states_idx + allpass_interpolator::n_states,
    impl_states_idx   = boxcar_states_idx + moving_average<2>::n_states,
    n_states          = impl_states_idx + Impl<1>::n_states
  };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (crange<vec_value_type_t<V>> c)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    static_assert (moving_average<2>::n_coeffs == 0, "Add initialization!");
    static_assert (Impl<1>::n_coeffs == 0, "Add initialization!");

    allpass_interpolator::init (c, vec_set<V> ((T) 0.5));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> s)
  {
    using T = vec_value_type_t<V>;

    // maybe memset delay line and boxcar?
    Impl<1>::template reset_states<V> (
      s.shrink_head (impl_states_idx * vec_traits<V>().size));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c,
    crange<vec_value_type_t<V>>       st,
    V                                 x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * n_states));

    // comments on the non-vectorized version
    V delayed = allpass_interpolator::tick (c, st, x);
    st        = st.shrink_head (allpass_interpolator::n_states * traits.size);

    V filtered = moving_average<2>::tick ({}, st, x);
    st         = st.shrink_head (moving_average<2>::n_states * traits.size);

    V eq = delayed + (delayed - filtered);
    return Impl<1>::tick ({}, st, eq);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
namespace detail {
template <uint order>
struct null_shaper {
  enum coeffs { n_coeffs };
  enum state { n_states };

  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c,
    crange<vec_value_type_t<V>>       st,
    V                                 x)
  {
    return x;
  }
};
} // namespace detail
//------------------------------------------------------------------------------
// initialization for the coefficients of all "fix_eq_and_delay" classes is the
// same, so "fix_eq_and_delay_coeff_initialization" can be used.
template <uint order>
struct fix_eq_and_delay_coeff_initialization
  : public fix_eq_and_delay<order, detail::null_shaper> {};

}} // namespace artv::adaa
