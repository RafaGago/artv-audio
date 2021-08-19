#pragma once

// Implementation heavily "inspired" on:
// https://github.com/jatinchowdhury18/ADAA

#include <cmath>

#include "artv-common/dsp/own/blocks/filters/moving_average.hpp"
#include "artv-common/dsp/own/blocks/misc/interpolators.hpp"

#include "artv-common/misc/short_ints.hpp"
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
  template <class T>
  static T tick (crange<const T>, crange<T>, T x)
  {
    return functions::fn (x);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
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
  template <class T>
  static void init_states (crange<T> s)
  {}
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_states_simd (crange<vec_value_type_t<V>> s)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    assert (st.size() >= n_states);

    T x_int   = functions::int_fn (x);
    T x1v     = st[x1];
    T x1v_int = st[x1_int];

    st[x1]     = x;
    st[x1_int] = x_int;

    T diff = x - x1v;
    if (abs (diff) >= epsilon (T {})) {
      return (x_int - x1v_int) / diff;
    }
    else {
      return (functions::fn ((x + x1v) * (T) 0.5));
    }
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V                           x)
  {
    // as this has to calculate both branches, it might not be worth bothering.
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (traits.size * n_states));

    T* x1v_ptr     = &st[x1 * traits.size];
    T* x1v_int_ptr = &st[x1_int * traits.size];

    V x1v     = vec_load<V> (x1v_ptr);
    V x1v_int = vec_load<V> (x1v_int_ptr);

    V x_int = functions::int_fn (x);

    vec_store (x1v_ptr, x);
    vec_store (x1v_int_ptr, x_int);

    V diff  = x - x1v;
    V big   = (x_int - x1v_int) / diff;
    V small = functions::fn ((x + x1v) * (T) 0.5);
    return (vec_abs (diff) >= epsilon (T {})) ? big : small;
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
  template <class T>
  static void init_states (crange<T> s)
  {}
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_states_simd (crange<vec_value_type_t<V>> s)
  {}
  //----------------------------------------------------------------------------
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
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
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
  template <class T>
  static void init (crange<T> c)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (c.size() >= n_coeffs);

    static_assert (moving_average<2>::n_coeffs == 0, "Add initialization!");
    static_assert (Impl<1>::n_coeffs == 0, "Add initialization!");

    allpass_interpolator::init (c, (T) 0.5);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_simd (crange<vec_value_type_t<V>> c)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    static_assert (moving_average<2>::n_coeffs == 0, "Add initialization!");
    static_assert (Impl<1>::n_coeffs == 0, "Add initialization!");

    allpass_interpolator::init_simd (c, vec_set<V> ((T) 0.5));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void init_states (crange<T> s)
  {
    static_assert (std::is_floating_point<T>::value, "");
    // maybe memset delay line and boxcar?
    Impl<1>::init_states (s.shrink_head (impl_states_idx));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_states_simd (crange<vec_value_type_t<V>> s)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");

    // maybe memset delay line and boxcar?
    Impl<1>::template init_states_simd<T, V> (
      s.shrink_head (impl_states_idx * vec_traits<V>().size));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T> c, crange<T> st, T x)
  {
    assert (c.size() >= n_coeffs);
    assert (st.size() >= n_states);
    // Explicit half sample delay to the main signal. The ADAA whaveshaper
    // of first order will add another half sample, as it behaves as a boxcar of
    // L=2. The total delay will be 1 sample.
    T delayed = allpass_interpolator::tick (c, st, x);
    st.shrink_head (allpass_interpolator::n_states);
    // Apply a boxcar to the input of the same size as the ADAA equivalent.
    // Delays half a sample
    T filtered = moving_average<2>::tick ({}, st, x);
    st.shrink_head (moving_average<2>::n_states);
    // The difference between the delayed main signal and the Boxcar output will
    // contain filtered away high frequencies. These are added as a pre-EQ, to
    // alleviate the ADAA filter HF loss with low or no oversampling.
    T eq = delayed + (delayed - filtered);
    return Impl<1>::tick ({}, st, eq);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c,
    crange<vec_value_type_t<V>>       st,
    V                                 x)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * n_states));

    // comments on the non-vectorized version
    V delayed = allpass_interpolator::tick_simd (c, st, x);
    st.shrink_head (allpass_interpolator::n_states * traits.size);

    V filtered = moving_average<2>::tick_simd ({}, st, x);
    st.shrink_head (moving_average<2>::n_states * traits.size);

    V eq = delayed + (delayed - filtered);
    return Impl<1>::tick_simd ({}, st, eq);
  }
  //----------------------------------------------------------------------------
};
namespace detail {
template <uint order>
struct null_shaper {
  enum coeffs { n_coeffs };
  enum state { n_states };

  template <class T>
  static T tick (crange<const T>, crange<T>, T)
  {}

  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c,
    crange<vec_value_type_t<V>>       st,
    V                                 x)
  {
    return x;
  }
};
} // namespace detail

// initialization for the coefficients of all "fix_eq_and_delay" classes is the
// same, so "fix_eq_and_delay_coeff_initialization" can be used.
template <uint order>
struct fix_eq_and_delay_coeff_initialization
  : public fix_eq_and_delay<order, detail::null_shaper> {};

}} // namespace artv::adaa
