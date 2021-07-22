#pragma once

// Implementation heavily "inspired" on:
// https://github.com/jatinchowdhury18/ADAA

#include <cmath>

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
  static simd_reg<T, simd_bytes> tick_aligned (
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
  static T tick (crange<const T>, crange<T> st, T x)
  {
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
  static simd_reg<T, simd_bytes> tick_aligned (
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
  static simd_reg<T, simd_bytes> tick_aligned (
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
}} // namespace artv::adaa
