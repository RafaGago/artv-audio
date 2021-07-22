#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

//------------------------------------------------------------------------------
namespace artv {

struct onepole_smoother {
  //----------------------------------------------------------------------------
  enum coeffs { b1, n_coeffs };
  enum state { z1, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void lowpass (crange<T> c, T freq, T sr)
  {
    static_assert (std::is_floating_point<T>::value, "");
    constexpr double pi_x2 = 6.283185307179586476925286766559;
    c[b1]                  = (T) exp (-pi_x2 * freq / sr);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       z, // state
    T               in)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (z.size() >= n_states);
    assert (c.size() >= n_coeffs);

    z[z1] = (in * (1. - c[b1])) + (z[z1] * c[b1]);
    return z[z1];
  }
  //----------------------------------------------------------------------------
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick (
    crange<const T>                                      c, // coeffs
    std::array<crange<T>, simd_reg<T, simd_bytes>::size> z, // state
    std::array<T, simd_reg<T, simd_bytes>::size>         in)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg             = simd_reg<T, simd_bytes>;
    constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= n_coeffs);
    for (uint i = 0; i < n_builtins; ++i) {
      assert (z[i].size() >= n_states);
      assert (z[i].size() >= n_states);
    }

    // I don't know if it's actually worth to use unaligned SIMD for this
    // actually.
    simdreg a0_v {1. - c[b1]};
    simdreg b1_v {c[b1]};
    simdreg in_v {in.data(), xsimd::unaligned_mode {}};
    simdreg z1_v;

    for (uint i = 0; i < n_builtins; ++i) {
      z1_v[i] = z[i][z1];
    }

    z1_v = (in_v * a0_v) + (z1_v * b1_v);

    for (uint i = 0; i < n_builtins; ++i) {
      z[i][z1] = z1_v[i];
    }
    return z1_v;
  }
  //----------------------------------------------------------------------------
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_aligned (
    crange<const T>         c, // coeffs, just 'b1'
    crange<T>               z, // state 'z1' 1 to N
    simd_reg<T, simd_bytes> in) // in' 1 to N
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg             = simd_reg<T, simd_bytes>;
    constexpr auto n_builtins = simdreg::size;

    assert (z.size() >= n_builtins * n_states);
    assert (c.size() >= n_coeffs);

    simdreg a0_v {((T) 1.) - c[b1]};
    simdreg b1_v {c[b1]};
    simdreg z1_v {z.data(), xsimd::aligned_mode {}};

    z1_v = (in * a0_v) + (z1_v * b1_v);

    z1_v.store_aligned (z.data());

    return z1_v;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// From ReEQ.
struct onepole {
  //----------------------------------------------------------------------------
  enum coeffs { b0, b1, a1, n_coeffs };
  enum state { z1, z0, n_states };
  //----------------------------------------------------------------------------
  static void lowpass (crange<double> co, double freq, double sr)
  {
    double w = tan (M_PI * freq / sr);
    double n = 1. / (1. + w);
    co[b0]   = w * n;
    co[b1]   = co[b0];
    co[a1]   = n * (w - 1.);
  }
  //----------------------------------------------------------------------------
  static void highpass (crange<double> co, double freq, double sr)
  {
    double w = tan (M_PI * freq / sr);
    double n = 1. / (1. + w);
    co[b0]   = n;
    co[b1]   = -co[b0];
    co[a1]   = n * (w - 1.);
  }
  //----------------------------------------------------------------------------
  static void repair_unsmoothable_coeffs (crange<double>, crange<const double>)
  {}
  //----------------------------------------------------------------------------
  static void reset_states (crange<double> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    st[z1] = in * co[b0] + st[z0] * co[b1] - st[z1] * co[a1];
    st[z0] = in;
    return st[z1];
  }
  //----------------------------------------------------------------------------
  static simd_dbl tick (
    crange<const double>          co, // coeffs
    std::array<crange<double>, 2> st, // state
    simd_dbl                      in)
  {
    assert (st.size() >= 2);
    assert (co.size() >= n_coeffs);
    assert (st[0].size() >= n_states);
    assert (st[1].size() >= n_states);

    // I don't know if it's worth to use SIMD for this actually.
    simd_dbl a1_v {co[a1]};
    simd_dbl b0_v {co[b0]};
    simd_dbl b1_v {co[b1]};
    simd_dbl z0_v {st[0][z0], st[1][z0]};
    simd_dbl z1_v {st[0][z1], st[1][z1]};

    z1_v      = (in * b0_v) + (z0_v * b1_v) - (z1_v * a1_v);
    st[0][z1] = z1_v[0];
    st[0][z0] = in[0];
    st[1][z1] = z1_v[1];
    st[1][z0] = in[1];
    return z1_v;
  }
};
//------------------------------------------------------------------------------
} // namespace artv
