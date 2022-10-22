#pragma once

#include <cmath>

#include "artv-common/dsp/own/parts/traits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"

namespace artv {
//------------------------------------------------------------------------------
// Not strictly a filter, just a linked low and high shelf. Move to another
// place?
class tilt_eq {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 2 * andy::svf::n_coeffs;
  static constexpr uint n_coeffs_int = 2 * andy::svf::n_coeffs_int;
  static constexpr uint n_states     = 2 * andy::svf::n_states;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl)
  {
    assert (co.size() >= n_coeffs);
    andy::svf::reset_coeffs (co, freq, q, gain_db, t_spl, lowshelf_tag {});
    co.cut_head (andy::svf::n_coeffs);
    andy::svf::reset_coeffs (co, freq, q, -gain_db, t_spl, highshelf_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    xspan<vec_value_type_t<V> const> co, // coeffs (single set)
    xspan<V>                         st, // states (interleaved, SIMD aligned)
    V                                in)
  {
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    in = andy::svf::tick (co, st, in);
    co.cut_head (andy::svf::n_coeffs);
    st.cut_head (andy::svf::n_states);
    return andy::svf::tick (co, st, in);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    xspan<V const> co, // states (interleaved, SIMD aligned)
    xspan<V>       st, // states (interleaved, SIMD aligned)
    V              in)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    in = andy::svf::tick (co, st, in);
    co.cut_head (andy::svf::n_coeffs);
    st.cut_head (andy::svf::n_states);
    return andy::svf::tick (co, st, in);
  }
  //----------------------------------------------------------------------------
};
// spectral tilt, from Faust's "spectral_tilt" sources. Tweaked to use bandwidth
// in octaves instead of absolute Hz. Useful for e.g get 3dB/oct filters. (alpha
// at 1/2). Requires very high orders to look decent on graphs.
template <uint N>
class spectral_tilt {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = N * (onepole_tdf2::n_coeffs + 1);
  static constexpr uint n_states     = N * onepole_tdf2::n_states;
  static constexpr uint n_coeffs_int = N * onepole_tdf2::n_coeffs_int;

  static_assert (N >= 2);
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co,
    V                   f0,
    V                   bw_oct,
    V                   alpha,
    vec_value_type_t<V> t_spl,
    raw_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V t     = vec_set<V> (t_spl);
    V w     = (T) 2 * (T) M_PI * f0;
    V f1    = f0 * vec_exp (bw_oct * (T) M_LN2);
    V srate = (T) 1 / t_spl;
    f1 = vec_min (f1, (srate * (T) 0.5) - (T) 10); // 10 Hz margin ont Nyquist
    V ratio = vec_pow (f1 / f0, vec_set<V> ((T) 1 / (T) (N - 1)));
    T c     = (T) 1 / tan ((T) 1 * (T) 0.5 * t_spl); // 1 = wd

    V* gains = &co[N * onepole_tdf2::n_coeffs];

    for (uint i = 0; i < N; ++i) {
      V mzh_v = mzh (w, ratio, alpha, t, i); // b0
      V mph_v = mph (w, ratio, t, i); // a0
      // 1 = b1
      gains[i] = mph / mzh;

      // BLT, (tf1s on faust). b1 = 1, b0 = mzh, a0 = mph, w1 = 1
      V d   = mph_v + c;
      V b1d = (mzh_v - (T) 1 * c) / d;
      V b0d = (mzh_v + (T) 1 * c) / d;
      V a1d = (mph_v - c) / d;

      onepole_tdf2::reset_coeffs (
        co.cut_head (onepole_tdf2::n_coeffs), a1d, b0d, b1d, raw_tag {});
    }
  }
  //----------------------------------------------------------------------------
  // "fc", the point of 0dB gain will be off unless high orders +20 are used.
  // Best effort.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co,
    V                   fc,
    V                   bw_oct,
    V                   alpha,
    vec_value_type_t<V> t_spl)
  {
    using T = vec_value_type_t<V>;
    // move back half number of octaves
    V f0 = fc * vec_exp (bw_oct * (T) -0.5 * (T) M_LN2);
    reset_coeffs (co, fc, bw_oct, alpha, t_spl, raw_tag {});

    V* gains = &co[N * onepole_tdf2::n_coeffs];

    // hack the last coefficient with extra gain offset to bring the response
    // from -gain to +gain.
    V half_db_range = bw_oct * -alpha * (6 + 6 / (T) (N - 1)) * (T) 0.5;
    gains[N - 1] *= vec_exp (half_db_range * (T) (1. / 20.) * (T) M_LN10);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (xspan<V const> co, xspan<V> st, V in)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V        out   = in;
    V const* gains = &co[N * onepole_tdf2::n_coeffs];

    for (uint i = 0; i < N; ++i) {
      auto c = co.cut_head (onepole_tdf2::n_coeffs);
      auto s = st.cut_head (onepole_tdf2::n_states);
      out    = gains[i] * onepole_tdf2::tick (c, s, out);
    }
    return out;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class V>
  static V prewarp (V w, V wp, V t)
  {
    using T = vec_value_type_t<V>;
    return wp * vec_tan (w * t * (T) 0.5) / vec_tan (wp * t * (T) 0.5);
  }
  //----------------------------------------------------------------------------
  template <class V>
  static V mz (V w, V ratio, V alpha, uint i)
  {
    using T = vec_value_type_t<V>;
    return w * vec_pow (ratio, (-alpha + (T) i));
  }
  //----------------------------------------------------------------------------
  template <class V>
  static V mp (V w, V ratio, uint i)
  {
    using T = vec_value_type_t<V>;
    return w * vec_pow (ratio, (T) i);
  }
  //----------------------------------------------------------------------------
  template <class V>
  static V mzh (V w, V ratio, V alpha, V t, uint i)
  {
    return prewarp (mz (w, ratio, alpha, i), w, t);
  }
  //----------------------------------------------------------------------------
  template <class V>
  static V mph (V w, V ratio, V t, uint i)
  {
    return prewarp (mp (w, ratio, i), w, t);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

} // namespace artv
