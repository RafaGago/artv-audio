#pragma once

#include <cassert>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/misc/simd.hpp"

//------------------------------------------------------------------------------
namespace artv {

// Classic IIR DC blocker based on a highpass
// https://ccrma.stanford.edu/~jos/filters/DC_Blocker.html
// https://www.musicdsp.org/en/latest/Filters/135-dc-filter.html?highlight=DC
struct iir_dc_blocker {
  //----------------------------------------------------------------------------
  enum coeffs { R, n_coeffs };
  enum state { x1, y1, n_states };
  //----------------------------------------------------------------------------
  // warning, if going to very low frequencies, use "double".
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c, V freq, vec_value_type_t<V> sr)
  {
    using T = vec_value_type_t<V>;
    assert (c.size() >= n_coeffs);

    c[R] = (T) 1. - ((T) M_PI * (T) 2. * freq / sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> c, crange<V> s, V x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= n_coeffs);
    return tick (s, x, c[R]);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const vec_value_type_t<V>> c, crange<V> s, V x)
  {
    assert (c.size() >= n_coeffs);
    return tick (s, x, vec_set<V> (c[R]));
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<V> s, V x, V Rv)
  {
    assert (s.size() >= n_states);

    auto y = x - s[x1] + Rv * s[y1];
    s[x1]  = x;
    s[y1]  = y;
    return y;
  }
};
//------------------------------------------------------------------------------
// I don't know if he is the original author but he does an extreme good job
// sharing knowledge and helping. Described here:
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=545280#top
struct mystran_dc_blocker {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs = onepole_smoother::n_coeffs;
  static constexpr uint n_states = onepole_smoother::n_states;
  //----------------------------------------------------------------------------
  // warning, if going to very low frequencies, use "double".
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c, V freq, vec_value_type_t<V> sr)
  {
    onepole_smoother::reset_coeffs (c, freq, sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> c, crange<V> s, V x)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    auto prev_lp_out = s[onepole_smoother::z1];
    onepole_smoother::tick (c, s, x);
    return x - prev_lp_out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c, // single coeff set
    crange<V>                         s,
    V                                 x)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    auto prev_lp_out = s[onepole_smoother::z1];
    onepole_smoother::tick (c, s, x);
    return x - prev_lp_out;
  }
  //----------------------------------------------------------------------------
};
// TODO: Things to try
// -Moving average DC blocker (linear phase):
//    https://www.dsprelated.com/showthread/comp.dsp/66509-1.php
//------------------------------------------------------------------------------
struct mystran_dc_blocker_2pole {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs = andy::svf_lowpass::n_coeffs;
  static constexpr uint n_states = andy::svf_lowpass::n_states + 1;
  //----------------------------------------------------------------------------
  // warning, if going to very low frequencies, use "double".
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c, V freq, vec_value_type_t<V> sr)
  {
    using T                   = vec_value_type_t<V>;
    constexpr T butterworth_q = M_SQRT1_2;

    andy::svf_lowpass::reset_coeffs (c, freq, vec_set<V> (butterworth_q), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> c, crange<V> s, V x)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    V  now         = andy::svf_lowpass::tick (c, s, x);
    V* prev        = s.shrink_head (andy::svf_lowpass::n_states).data();
    V  prev_lp_out = *prev;
    *prev          = now;
    return x - prev_lp_out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c, // single coeff set
    crange<V>                         s,
    V                                 x)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    V  now         = andy::svf_lowpass::tick (c, s, x);
    V* prev        = s.shrink_head (andy::svf_lowpass::n_states).data();
    V  prev_lp_out = *prev;
    *prev          = now;
    return x - prev_lp_out;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
