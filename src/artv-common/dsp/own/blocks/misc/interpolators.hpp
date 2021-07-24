#pragma once

#include "artv-common/misc/util.hpp"

namespace artv {

// Seen on Saike TanhAA saturator
struct allpass_interpolator {
  //----------------------------------------------------------------------------
  enum coeffs { nu, n_coeffs };
  enum state { y0, d0, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void init (crange<T> c, T frac)
  {
    c[nu] = (1.0 - frac) / (1.0 + frac);
  }
  //----------------------------------------------------------------------------
  template <uint simd_bytes, class T>
  static void init_multi_aligned (crange<T> c, simd_reg<T, simd_bytes> frac)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= (n_coeffs * n_builtins));

    auto v = (simdreg {(T) 1.0} - frac) / (1.0 + frac);
    v.store_aligned (&c[nu * n_builtins]);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr T tick (
    crange<const T> c, // coeffs
    crange<T>       z, // state
    T               in)
  {
    z[y0] = c[nu] * in + z[d0] - c[nu] * z[y0];
    z[d0] = in;
    return z[y0];
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>         c, // coeffs interleaved, ready to SIMD load
    crange<T>               z, // states interleaved, ready to SIMD load
    simd_reg<T, simd_bytes> in) // N inputs ready to SIMD load
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= n_builtins * n_coeffs);
    assert (z.size() >= n_builtins * n_states);

    T const* nu_ptr = &c[nu * n_builtins];
    T*       d0_ptr = &z[d0 * n_builtins];
    T*       y0_ptr = &z[y0 * n_builtins];

    simdreg nu_v {nu_ptr, xsimd::aligned_mode {}};
    simdreg d0_v {d0_ptr, xsimd::aligned_mode {}};
    simdreg y0_v {y0_ptr, xsimd::aligned_mode {}};

    y0_v = nu_v * in + d0_v - nu_v * y0_v;

    y0_v.store_aligned (y0_ptr);
    in.store_aligned (d0_ptr);
    return y0_v;
  }
};
//----------------------------------------------------------------------------

} // namespace artv
