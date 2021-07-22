#pragma once

#include <cassert>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// a slew limiter. Notice that an envelope follower can be made by modifying
// the input, and sending:
//
// tick (fabs (in)): peak follower
// tick (in * in): mean abs
// sqrt (tick (in * in)): instant RMS
//
struct slew_limiter {
  //----------------------------------------------------------------------------
  enum coeffs { attack, release, n_coeffs };
  enum state { prev, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void init (crange<T> c, T attack_sec, T release_sec, T samplerate)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (attack_sec >= 0.);
    assert (release_sec >= 0.);

#if 0 // TODO: xsimd broken with -ffast-math
    simd_dbl v {attack_sec, release_sec};
    v *= samplerate;
    v          = xsimd::exp (simd_dbl {-1.} / v);
    c[attack]  = v[0];
    c[release] = v[1];
#else
    c[attack]  = exp (-1. / (samplerate * attack_sec));
    c[release] = exp (-1. / (samplerate * release_sec));
#endif
  }
  //----------------------------------------------------------------------------
  template <size_t simd_bytes, class T>
  static void init_multi_aligned (
    crange<T>               c,
    simd_reg<T, simd_bytes> attack_sec,
    simd_reg<T, simd_bytes> release_sec,
    T                       samplerate)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= (n_coeffs * n_builtins));

#if 0 // TODO: xsimd broken with -ffast-math
    attack_sec *= samplerate;
    attack_sec = xsimd::exp (simd_dbl {1.} / attack_sec);
    attack_sec.store_aligned (&c[attack * n_builtins]);

    release_sec *= samplerate ;
    release_sec = xsimd::exp (simd_dbl {1.} / release_sec);
    release_sec.store_aligned (&c[release * n_builtins]);
#else
    for (uint i = 0; i < n_builtins; ++i) {
      c[(attack * n_builtins) + i]  = exp (-1. / (samplerate * attack_sec[i]));
      c[(release * n_builtins) + i] = exp (-1. / (samplerate * release_sec[i]));
    }
#endif
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               in)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (s.size() >= n_states);
    assert (c.size() >= n_coeffs);

    T in_prev = s[prev];
    T coeff   = in_prev < in ? c[attack] : c[release];
    s[prev]   = in + coeff * (in_prev - in);
    return s[prev];
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>         c, // coeffs interleaved, ready to SIMD load
    crange<T>               s, // states interleaved, ready to SIMD load
    simd_reg<T, simd_bytes> in)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= n_builtins * n_coeffs);
    assert (s.size() >= n_builtins * n_states);

    simdreg in_prev, att, rel;
    in_prev.load_aligned (&s[prev * n_builtins]);
    att.load_aligned (&c[attack * n_builtins]);
    rel.load_aligned (&c[release * n_builtins]);

    auto coeffs = xsimd::select (in_prev < in, att, rel);
    in_prev     = in + coeffs * (in_prev - in);
    in_prev.store_aligned (&s[prev * n_builtins]);
    return in_prev;
  }
};
//------------------------------------------------------------------------------

} // namespace artv
