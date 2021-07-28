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
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>,
    crange<T>,
    simd_reg<T, simd_bytes> x)
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
  template <size_t simd_bytes, class T>
  static void init_states_multi_aligned (crange<T> s)
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
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>,
    crange<T>               st,
    simd_reg<T, simd_bytes> x)
  {
    // as this has to calculate both branches, it might not be worth bothering.
    static_assert (std::is_floating_point<T>::value, "");
    using regtype                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = regtype::size;

    assert (st.size() >= n_builtins * n_states);

    T* x1v_ptr     = &st[x1 * n_builtins];
    T* x1v_int_ptr = &st[x1_int * n_builtins];

    regtype x1v, x1v_int;
    x1v.load_aligned (x1v_ptr);
    x1v_int.load_aligned (x1v_int_ptr);

    regtype x_int = functions::int_fn (x);

    x.store_aligned (x1v_ptr);
    x_int.store_aligned (x1v_int_ptr);

    regtype diff  = x - x1v;
    regtype big   = (x_int - x1v_int) / diff;
    regtype small = functions::fn ((x + x1v) * (T) 0.5);
    return xsimd::select (
      xsimd::abs (diff) >= regtype {epsilon (T {})}, big, small);
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
  template <size_t simd_bytes, class T>
  static void init_states_multi_aligned (crange<T> s)
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
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>,
    crange<T>               st,
    simd_reg<T, simd_bytes> x)
  {
    // Real SIMD: TODO. Might not be worth because of the high number of
    // branches.
    static_assert (!std::is_same_v<T, T>, "To be implemented");
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
  template <size_t simd_bytes, class T>
  static void init_multi_aligned (crange<T> c)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= (n_coeffs * n_builtins));

    static_assert (moving_average<2>::n_coeffs == 0, "Add initialization!");
    static_assert (Impl<1>::n_coeffs == 0, "Add initialization!");

    allpass_interpolator::init_multi_aligned<simd_bytes, T> (
      c, simdreg {(T) 0.5});
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
  template <size_t simd_bytes, class T>
  static void init_states_multi_aligned (crange<T> s)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;
    // maybe memset delay line and boxcar?
    Impl<1>::template init_states_multi_aligned<simd_bytes, T> (
      s.shrink_head (impl_states_idx * n_builtins));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T> c, crange<T> st, T x)
  {
    assert (c.size() >= n_coeffs);
    assert (st.size() >= n_states);

    // Explicit half sample delay to the main signal, the ADAA whaveshaper
    // of first order will add another half sample, as it behaves as a boxcar of
    // L=2. This makes the delay to happen at sample boundaries.
    T delayed = allpass_interpolator::tick (c, st, x);
    st.shrink_head (allpass_interpolator::n_states);
    // Apply boxcar to the input, delays half sample and rolls of highs
    T filtered = moving_average<2>::tick ({}, st, x);
    st.shrink_head (moving_average<2>::n_states);
    // As the ADAA implementation behaves as a boxcar, the input signal is
    // pre-EQ'd to result in a flat frequency response after the ADAA process
    // runs.
    T eq = delayed + (delayed - filtered);
    return Impl<1>::tick ({}, st, eq);
  }
  //----------------------------------------------------------------------------
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>         c,
    crange<T>               st,
    simd_reg<T, simd_bytes> x)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using regtype                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = regtype::size;

    assert (c.size() >= (n_builtins * n_coeffs));
    assert (st.size() >= (n_builtins * n_states));

    // comments on the non-vectorized version
    regtype delayed
      = allpass_interpolator::tick_multi_aligned<simd_bytes, T> (c, st, x);
    st.shrink_head (allpass_interpolator::n_states * n_builtins);

    regtype filtered
      = moving_average<2>::tick_multi_aligned<simd_bytes, T> ({}, st, x);
    st.shrink_head (moving_average<2>::n_states * n_builtins);

    regtype eq = delayed + (delayed - filtered);
    return Impl<1>::template tick_multi_aligned<simd_bytes, T> ({}, st, eq);
  }
  //----------------------------------------------------------------------------
};
namespace detail {
template <uint order>
struct null_shaper {
  enum coeffs { n_coeffs };
  enum state { n_states };

  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {}

  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>,
    crange<T>,
    simd_reg<T, simd_bytes>)
  {}
};
} // namespace detail

// initialization for the coefficients of all "fix_eq_and_delay" classes is the
// same, so "fix_eq_and_delay_coeff_initialization" can be used.
template <uint order>
struct fix_eq_and_delay_coeff_initialization
  : public fix_eq_and_delay<order, detail::null_shaper> {};

}} // namespace artv::adaa
