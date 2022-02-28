#pragma once

// Implementation heavily "inspired" on:
// https://github.com/jatinchowdhury18/ADAA

#include <cmath>

#include "artv-common/dsp/own/parts/filters/moving_average.hpp"
#include "artv-common/dsp/own/parts/misc/interpolators.hpp"

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv { namespace adaa {
//------------------------------------------------------------------------------
template <class T>
static constexpr float epsilon (T)
{
  // untested
  return (T) 1.0e-5;
}
//------------------------------------------------------------------------------
template <class functions, uint order>
class waveshaper;
//------------------------------------------------------------------------------
template <class functions>
class waveshaper<functions, 0> {
public:
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, V x)
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
  enum coeffs_int { n_coeffs_int };
  enum state { x1, x1_int, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V> st, V x)
  {
    // as this has to calculate both branches, it might not be worth bothering.
    using T = vec_value_type_t<V>;

    assert (st.size() >= n_states);

    V x1v     = st[x1];
    V x1v_int = st[x1_int];

    V x_int = functions::int_fn (x);

    st[x1]     = x;
    st[x1_int] = x_int;

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
  enum coeffs_int { n_coeffs_int };
  enum state { x1, x2, x2_der, x1_int2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
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
  static V tick (crange<const V>, crange<V> st, V x)
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
template <uint order, template <uint, class...> class Impl, class... Ts>
class fix_eq_and_delay;

template <template <uint, class...> class Impl, class... Ts>
class fix_eq_and_delay<1, Impl, Ts...> {
public:
  enum coeffs {
    allpass_coeffs_idx,
    boxcar_coeffs_idx = allpass_coeffs_idx + allpass_interpolator::n_coeffs,
    impl_coeffs_idx   = boxcar_coeffs_idx + moving_average<2>::n_coeffs,
    n_coeffs          = impl_coeffs_idx + Impl<1, Ts...>::n_coeffs
  };
  enum coeffs_int { n_coeffs_int };
  enum state {
    allpass_states_idx,
    boxcar_states_idx = allpass_states_idx + allpass_interpolator::n_states,
    impl_states_idx   = boxcar_states_idx + moving_average<2>::n_states,
    n_states          = impl_states_idx + Impl<1, Ts...>::n_states
  };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c)
  {
    using T = vec_value_type_t<V>;

    assert (c.size() >= n_coeffs);

    static_assert (moving_average<2>::n_coeffs == 0, "Add initialization!");
    static_assert (Impl<1, Ts...>::n_coeffs == 0, "Add initialization!");

    allpass_interpolator::reset_coeffs (c, vec_set<V> ((T) 0.5));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> s)
  {
    // maybe memset delay line and boxcar?
    Impl<1, Ts...>::template reset_states<V> (s.advanced (impl_states_idx));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> c, crange<V> st, V x)
  {
    assert (c.size() >= n_coeffs);
    assert (st.size() >= n_states);

    // comments on the non-vectorized version
    V delayed = allpass_interpolator::tick (c, st, x);
    st.cut_head (allpass_interpolator::n_states);

    V filtered = moving_average<2>::tick ({}, st, x);
    st.cut_head (moving_average<2>::n_states);

    V eq = delayed + (delayed - filtered);
    return Impl<1, Ts...>::tick ({}, st, eq);
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
  static V tick (crange<const V> c, crange<V> st, V x)
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
