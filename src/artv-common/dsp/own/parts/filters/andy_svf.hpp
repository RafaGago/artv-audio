#pragma once

/* Fiters from: https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf*/

#include <array>
#include <cmath>
#include <gcem.hpp>
#include <type_traits>

#include "artv-common/dsp/own/parts/filters/zdf.hpp"
#include "artv-common/dsp/own/parts/traits.hpp"

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace andy {
namespace detail {
//------------------------------------------------------------------------------
template <class T>
struct svf_coeffs {
  T a1;
  T a2;
  T a3;
  T k;
};

template <class V>
static svf_coeffs<V> get_main_coeffs (V g, V q)
{
  // "T" should always be a builtin type.
  using T = vec_value_type_t<V>;

  svf_coeffs<V> ret;
  ret.k  = (T) 1.0 / q;
  ret.a1 = (T) 1.0 / ((T) 1.0 + g * (g + ret.k));
  ret.a2 = g * ret.a1;
  ret.a3 = g * ret.a2;
  return ret;
}

template <class V, class T>
static svf_coeffs<V> get_main_coeffs (V freq, V q, T t_spl)
{
  // "T" should always be a builtin type.
  return get_main_coeffs (vec_tan ((T) M_PI * freq * t_spl), q);
}
//------------------------------------------------------------------------------
template <class V>
struct svf_coeffs_ext : public svf_coeffs<V> {
  V A;
};

static constexpr uint bell_flag   = 1 << 0;
static constexpr uint lshelf_flag = 1 << 1;
static constexpr uint hshelf_flag = 1 << 2;

template <class V, class T>
static svf_coeffs_ext<V> get_main_coeffs (
  V    freq,
  V    q,
  V    db,
  T    t_spl,
  uint flags = 0)
{
  // "T" should always be a builtin type.
  static_assert (std::is_floating_point<T>::value, "");

  svf_coeffs_ext<V> ret;
  V                 g;

  if constexpr (is_vec_v<V>) {
    ret.A = vec_exp (db * (T) (1. / 40.) * (T) M_LN10);
    g     = vec_tan ((T) M_PI * freq * t_spl);
  }
  else {
    ret.A = exp (db * (T) (1. / 40.) * (T) M_LN10);
    g     = tan ((T) M_PI * freq * t_spl);
  }
  if ((flags & bell_flag)) {
    q *= ret.A;
  }
  if ((flags & lshelf_flag)) {
    g /= vec_sqrt (ret.A);
  }
  if ((flags & hshelf_flag)) {
    g *= vec_sqrt (ret.A);
  }
  ret.k  = (T) 1.0 / q;
  ret.a1 = (T) 1.0 / ((T) 1.0 + g * (g + ret.k));
  ret.a2 = g * ret.a1;
  ret.a3 = g * ret.a2;
  return ret;
}

// https://cytomic.com/files/dsp/SvfLinearTrapezoidalSin.pdf
// UNTESTED!
template <class V, class T>
static svf_coeffs_ext<V> get_main_coeffs_precise (
  V    freq,
  V    q,
  V    db,
  T    t_spl,
  uint flags = 0)
{
  // "T" should always be a builtin type.
  static_assert (std::is_floating_point<T>::value, "");

  svf_coeffs_ext<V> ret;
  V                 s1, s2;

  V w = (T) M_PI * freq * t_spl;
  if constexpr (is_vec_v<V>) {
    ret.A = vec_exp (db * (T) (1. / 40.) * (T) M_LN10);
    s1    = vec_sin (w);
    // sometimes sin and cos can be returned on the same op.
    // sin(2x) = 2*cos(x)*sin(x)
    s2 = (T) 2 * vec_cos (w) * s1;
  }
  else {
    ret.A = exp (db * (T) (1. / 40.) * (T) M_LN10);
    s1    = vec_sin (w);
    // sometimes sin and cos can be returned on the same op.
    // sin(2x) = 2*cos(x)*sin(x)
    s2 = (T) 2 * vec_cos (w) * s1;
  }
  V nrm = (T) 1 / ((T) 2 * s2);

  ret.k  = (T) 1 / q;
  ret.a1 = s2 * nrm;
  ret.a2 = ((T) -2 * s1 * s1 - ret.k * s2) * nrm;
  ret.a3 = ((T) 2 * s1 * s1) * nrm;
  return ret;
}
//------------------------------------------------------------------------------
template <class T>
struct svf_tick_result {
  T v1;
  T v2;
  T s1;
  T s2;
};

template <class T>
static svf_tick_result<T> svf_tick (T s1, T s2, T a1, T a2, T a3, T v0)
{
  svf_tick_result<T> ret;

  T v3   = v0 - s2;
  ret.v1 = a1 * s1 + a2 * v3;
  ret.v2 = s2 + a2 * s1 + a3 * v3;
  ret.s1 = (ret.v1 * 2.) - s1;
  ret.s2 = (ret.v2 * 2.) - s2;
  return ret;
}

//------------------------------------------------------------------------------
// Enables the SVF multimode output. It is compile-time optimized depending on
// the enabled modes. Use any of the "_tag" classes as template parameters to
// enable a given mode. The outputs are returned on the "tick" function on the
// order on which the tags appear (natural order, left to right). If only one
// tag is passed the output is returned as a single value instead of an array.
//------------------------------------------------------------------------------
template <class Derived, class... Tags>
struct svf_multimode {

  // maybe do tag validation, (things won't compile anyways).
  using enabled_modes = mp11::mp_unique<mp11::mp_list<Tags...>>;

  template <class tag>
  static constexpr bool contains = mp11::mp_contains<enabled_modes, tag>::value;

  // clang-format off
  static constexpr bool needs_k_coeff =
    contains<highpass_tag> ||
    contains<notch_tag> ||
    contains<peak_tag> ||
    contains<bandpass_tag> ||
    contains<allpass_tag>;
  // clang-format on

  static constexpr bool returns_array = mp11::mp_size<enabled_modes>::value > 1;

  static_assert (
    mp11::mp_size<enabled_modes>::value >= 1,
    "One mode needs to be enabled");
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, k, n_coeffs = k + (uint) needs_k_coeff };
  enum coeffs_int { n_coeffs_int };
  enum state { s1, s2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> c, V freq, V q, vec_value_type_t<V> t_spl)
  {
    assert (c.size() >= n_coeffs);

    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    if constexpr (needs_k_coeff) {
      c[k] = coeffs.k;
    }
    if constexpr (!std::is_same_v<Derived, void>) {
      Derived::template after_reset_coeffs (c, coeffs);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    no_prewarp)
  {
    assert (c.size() >= n_coeffs);
    using T = vec_value_type_t<V>;

    auto coeffs = detail::get_main_coeffs ((T) M_PI * freq * t_spl, q);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    if constexpr (needs_k_coeff) {
      c[k] = coeffs.k;
    }
    if constexpr (!std::is_same_v<Derived, void>) {
      Derived::template after_reset_coeffs (c, coeffs);
    }
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
  static auto tick (xspan<const V> co, xspan<V> st, V in)
  {
    return tick_impl<V, V> (co, st, in);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (xspan<const vec_value_type_t<V>> co, xspan<V> st, V in)
  {
    return tick_impl<V, vec_value_type_t<V>> (co, st, in);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class V, class VT, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick_impl (
    xspan<const VT> co, // coeffs (V builtin type (single set) or V (SIMD))
    xspan<V>        st, // states (interleaved, SIMD aligned)
    V               v0)
  {
    using T = vec_value_type_t<V>;
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    V a1_, a2_, a3_, k_;
    if constexpr (std::is_same_v<std::remove_const_t<VT>, V>) {
      a1_ = co[a1];
      a2_ = co[a2];
      a3_ = co[a3];
      if constexpr (needs_k_coeff) {
        k_ = co[k];
      }
    }
    else if constexpr (std::is_same_v<std::remove_const_t<VT>, T>) {
      a1_ = vec_set<V> (co[a1]);
      a2_ = vec_set<V> (co[a2]);
      a3_ = vec_set<V> (co[a3]);
      if constexpr (needs_k_coeff) {
        k_ = vec_set<V> (co[k]);
      }
    }
    else {
      static_assert (
        sizeof (V) == 0, "T type must be either V or V's builtin type");
    }

    auto tick_r = detail::svf_tick (st[s1], st[s2], a1_, a2_, a3_, v0);
    st[s1]      = tick_r.s1;
    st[s2]      = tick_r.s2;

    std::array<V, mp11::mp_size<enabled_modes>::value> ret;

    mp_foreach_idx (enabled_modes {}, [&] (auto index, auto mode) {
      using tag = decltype (mode);

      if constexpr (std::is_same_v<tag, lowpass_tag>) {
        ret[index.value] = tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, highpass_tag>) {
        ret[index.value] = v0 - k_ * tick_r.v1 - tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, bandpass_q_gain_tag>) {
        ret[index.value] = tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, bandpass_tag>) {
        ret[index.value] = tick_r.v1 * k_;
      }
      else if constexpr (std::is_same_v<tag, notch_tag>) {
        ret[index.value] = v0 - k_ * tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, peak_tag>) {
        ret[index.value] = v0 - k_ * tick_r.v1 - (T) 2. * tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, allpass_tag>) {
        ret[index.value] = v0 - (T) 2. * k_ * tick_r.v1;
      }
    });

    if constexpr (returns_array) {
      return ret;
    }
    else {
      return ret[0];
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
template <class... Tags>
using svf_multimode = detail::svf_multimode<void, Tags...>;

using svf_lowpass  = svf_multimode<lowpass_tag>;
using svf_highpass = svf_multimode<highpass_tag>;
using svf_allpass  = svf_multimode<allpass_tag>;
using svf_bandpass = svf_multimode<bandpass_tag>;
using svf_notch    = svf_multimode<notch_tag>;
using svf_peak     = svf_multimode<peak_tag>;
//------------------------------------------------------------------------------
// 1 pole adding the ability to get the S and G parameters for zdf feedback
// calculations. Using a separate class because it needs some extra coefficient
// memory to avoid divisions
template <class... Tags>
class svf_multimode_zdf
  : private detail::svf_multimode<svf_multimode_zdf<Tags...>, Tags...> {
public:
  using base = detail::svf_multimode<svf_multimode_zdf<Tags...>, Tags...>;

  template <class Tag>
  static constexpr bool contains = base::template contains<Tag>;

  // clang-format off
  static constexpr bool needs_b1_coeff =
    contains<highpass_tag> ||
    contains<bandpass_tag> ||
    contains<notch_tag> ||
    contains<peak_tag> ||
    contains<allpass_tag>;

  static constexpr bool needs_b2_coeff =
    contains<lowpass_tag> ||
    contains<bandpass_tag> ||
    contains<notch_tag> ||
    contains<peak_tag> ||
    contains<allpass_tag>;
  // clang-format on

  enum coeffs {
    b1       = base::n_coeffs,
    b2       = base::n_coeffs + (int) needs_b1_coeff,
    n_coeffs = b2 + (int) needs_b2_coeff
  };
  enum coeffs_int { n_coeffs_int = base::n_coeffs_int };
  enum state { n_states = base::n_states };

  static constexpr uint n_modes_enabled
    = mp11::mp_size<typename base::enabled_modes>::value;
  // coeffs are G and S pairs
  static constexpr uint n_zdf_coeffs = n_modes_enabled * 2;
  //----------------------------------------------------------------------------
  using base::reset_coeffs;
  using base::reset_states;
  using base::tick;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (xspan<const V> co, xspan<V> st, zdf::gs_coeffs_tag)
  {
    std::array<V, n_zdf_coeffs> ret;
    tick_impl<V, V> (co, st, ret, zdf::gs_coeffs_tag {});
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    xspan<const vec_value_type_t<V>> co,
    xspan<V>                         st,
    zdf::gs_coeffs_tag)
  {
    std::array<V, n_zdf_coeffs> ret;
    tick_impl<V, vec_value_type_t<V>> (co, st, ret, zdf::gs_coeffs_tag {});
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    xspan<const V> co,
    xspan<V>       st,
    xspan<V>       G_S,
    zdf::gs_coeffs_tag)
  {
    tick_impl<V, V> (co, st, G_S, zdf::gs_coeffs_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    xspan<const vec_value_type_t<V>> co,
    xspan<V>                         st,
    xspan<V>                         G_S,
    zdf::gs_coeffs_tag)
  {
    tick_impl<V, vec_value_type_t<V>> (co, st, G_S, zdf::gs_coeffs_tag {});
  }
  //----------------------------------------------------------------------------
private:
  // return G and S for each mode
  template <class V, class VT, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick_impl (
    xspan<const VT> co, // coeffs (V builtin type (single set) or V (SIMD))
    xspan<const V>  st, // states (interleaved, SIMD aligned)
    xspan<V>        G_S,
    zdf::gs_coeffs_tag)
  {
    using T = vec_value_type_t<V>;
    assert (G_S.size() >= n_zdf_coeffs);
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V a1_, a2_, a3_, b1_, b2_;
    if constexpr (std::is_same_v<VT, V>) {
      a1_ = co[base::a1];
      a2_ = co[base::a2];
      a3_ = co[base::a3];
      if constexpr (needs_b1_coeff) {
        b1_ = co[b1];
      }
      if constexpr (needs_b2_coeff) {
        b2_ = co[b2];
      }
    }
    else if constexpr (std::is_same_v<VT, T>) {
      a1_ = vec_set<V> (co[base::a1]);
      a2_ = vec_set<V> (co[base::a2]);
      a3_ = vec_set<V> (co[base::a3]);
      if constexpr (needs_b1_coeff) {
        b1_ = vec_set<V> (co[b1]);
      }
      if constexpr (needs_b2_coeff) {
        b2_ = vec_set<V> (co[b2]);
      }
    }
    else {
      static_assert (
        sizeof (V) == 0, "T type must be either V or V's builtin type");
    }

    mp_foreach_idx (
      typename base::enabled_modes {}, [&] (auto index, auto mode) {
        using tag    = decltype (mode);
        uint const i = index.value * 2;

        // b1 = k * a1, b2 = gk * a1 = k * a2
        // formulas below computed with
        // misc/sympy/2order-tpt.py

        if constexpr (std::is_same_v<tag, lowpass_tag>) {
          //         2
          //        g ⋅x + g⋅k⋅s₂ + g⋅s₁ + s₂
          // LP_y = ─────────────────────────
          //                2
          //               g  + g⋅k + 1
          //              2
          //             g
          // LP_G = ────────────
          //         2
          //        g  + g⋅k + 1
          //        g⋅k⋅s₂ + g⋅s₁ + s₂
          // LP_S = ──────────────────
          //            2
          //           g  + g⋅k + 1

          G_S[i + 0] = a3_; // G
          G_S[i + 1] = st[base::s1] * a2_ + st[base::s2] * (a1_ + b2_); // S
        }
        else if constexpr (std::is_same_v<tag, highpass_tag>) {
          //        -s₁⋅(g + k) - s₂ + x
          // HP_y = ────────────────────
          //             2
          //            g  + g⋅k + 1
          //             1
          // HP_G = ────────────
          //         2
          //        g  + g⋅k + 1
          //        -(g⋅s₁ + k⋅s₁ + s₂)
          // HP_S = ────────────────────
          //             2
          //            g  + g⋅k + 1

          G_S[i + 0] = a1_; // G
          G_S[i + 1] = st[base::s1] * (a2_ + b1_) + st[base::s2] * a1_;
          G_S[i + 1] = -G_S[i + 1]; // S
        }
        else if constexpr (std::is_same_v<tag, bandpass_q_gain_tag>) {
          //            -g⋅(s₂ - x) + s₁
          // BP_raw_y = ────────────────
          //               2
          //              g  + g⋅k + 1
          //                 g
          // BP_raw_G = ────────────
          //             2
          //            g  + g⋅k + 1
          //             -g⋅s₂ + s₁
          // BP_raw_S = ────────────
          //             2
          //            g  + g⋅k + 1
          G_S[i + 0] = a2_; // G
          G_S[i + 1] = st[base::s1] * a1_ - st[base::s2] * a2_; // S
        }
        else if constexpr (std::is_same_v<tag, bandpass_tag>) {
          //        -k⋅(g⋅(s₂ - x) - s₁)
          // BP_y = ─────────────────────
          //              2
          //             g  + g⋅k + 1
          //            g⋅k
          // BP_G = ────────────
          //         2
          //        g  + g⋅k + 1
          //        k⋅(-g⋅s₂ + s₁)
          // BP_S = ──────────────
          //          2
          //         g  + g⋅k + 1
          G_S[i + 0] = b2_; // G
          G_S[i + 1] = st[base::s1] * b1_ - st[base::s2] * b2_; // S
        }
        else if constexpr (std::is_same_v<tag, notch_tag>) {
          //            2
          //           g ⋅x + g⋅k⋅s₂ - k⋅s₁ + x
          // Notch_y = ────────────────────────
          //                  2
          //                 g  + g⋅k + 1
          //               2
          //              g  + 1
          // Notch_G = ────────────
          //            2
          //           g  + g⋅k + 1
          //           k⋅(g⋅s₂ - s₁)
          // Notch_S = ─────────────
          //             2
          //            g  + g⋅k + 1
          G_S[i + 0] = a3_ + a1_; // G
          G_S[i + 1] = st[base::s2] * b2_ - st[base::s1] * b1_; // S
        }
        else if constexpr (std::is_same_v<tag, peak_tag>) {
          //           2
          //          g ⋅x + g⋅k⋅s₂ + 2⋅g⋅s₁ + k⋅s₁ + 2⋅s₂ - x
          // Peak_y = ────────────────────────────────────────
          //                         2
          //                        g  + g⋅k + 1
          //              2
          //             g  - 1
          // Peak_G = ────────────
          //           2
          //          g  + g⋅k + 1
          //          g⋅k⋅s₂ + 2⋅g⋅s₁ + k⋅s₁ + 2⋅s₂
          // Peak_S = ─────────────────────────────
          //                    2
          //                   g  + g⋅k + 1
          G_S[i + 0] = a3_ - a1_; // G
          G_S[i + 1] = st[base::s1] * ((T) 2 * a2_ + b1_);
          G_S[i + 1] += st[base::s2] * ((T) 2 * a1_ + b2_);
        }
        else if constexpr (std::is_same_v<tag, allpass_tag>) {
          //                                 ⎛ 2          ⎞
          //        2⋅k⋅(g⋅(s₂ - x) - s₁) + x⋅⎝g  + g⋅k + 1⎠
          // AP_y = ────────────────────────────────────────
          //                       2
          //                      g  + g⋅k + 1
          //         2
          //        g  - g⋅k + 1
          // AP_G = ────────────
          //         2
          //        g  + g⋅k + 1
          //        2⋅k⋅(g⋅s₂ - s₁)
          // AP_S = ───────────────
          //           2
          //          g  + g⋅k + 1
          G_S[i + 0] = a3_ - b2_ + a1_; // G
          G_S[i + 1] = st[base::s2] * b2_ - st[base::s1] * b1_;
          G_S[i + 1] *= (T) 2; // S
        }
      });
  }
  //----------------------------------------------------------------------------
  friend class detail::svf_multimode<svf_multimode_zdf<Tags...>, Tags...>;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static inline void after_reset_coeffs (xspan<V> co, detail::svf_coeffs<V> c)
  {
    assert (co.size() >= n_coeffs);
    if constexpr (needs_b1_coeff) {
      //             k
      // b1 = ──────────────
      //        2
      //       g  + g⋅k + 1
      co[b1] = c.a1 * c.k;
    }

    if constexpr (needs_b2_coeff) {
      //         g * k
      // b2 = ──────────────
      //       2
      //      g  + g⋅k + 1
      co[b2] = c.a2 * c.k;
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
using svf_zdf_lowpass  = svf_multimode_zdf<lowpass_tag>;
using svf_zdf_highpass = svf_multimode_zdf<highpass_tag>;
using svf_zdf_allpass  = svf_multimode_zdf<allpass_tag>;
using svf_zdf_bandpass = svf_multimode_zdf<bandpass_tag>;
using svf_zdf_notch    = svf_multimode_zdf<notch_tag>;
using svf_zdf_peak     = svf_multimode_zdf<peak_tag>;
//------------------------------------------------------------------------------
// less eficient version that can switch modes on the fly. Useful for e.g. EQs.
struct svf {
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, m0, m1, m2, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { s1, s2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    lowpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = vec_set<V> ((T) 0.);
    c[m2] = vec_set<V> ((T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    highpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = -coeffs.k;
    c[m2] = vec_set<V> ((T) -1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    bandpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    assert (c.size() >= n_coeffs);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = coeffs.k;
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    bandpass_q_gain_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    assert (c.size() >= n_coeffs);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = vec_set<V> ((T) 1.);
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    peak_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    assert (c.size() >= n_coeffs);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = -coeffs.k;
    c[m2] = vec_set<V> ((T) -2.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    allpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = (T) -2. * coeffs.k;
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    notch_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, t_spl);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = (T) coeffs.k;
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> t_spl,
    bell_tag)
  {
    using T = vec_value_type_t<V>;
    auto coeffs
      = detail::get_main_coeffs (freq, q, db, t_spl, detail::bell_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = coeffs.k * (coeffs.A * coeffs.A - (T) 1.);
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> t_spl,
    bell_bandpass_tag)
  {
    using T = vec_value_type_t<V>;
    auto coeffs
      = detail::get_main_coeffs (freq, q, db, t_spl, detail::bell_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = coeffs.k * (coeffs.A * coeffs.A - (T) 1.);
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> t_spl,
    lowshelf_tag)
  {
    using T = vec_value_type_t<V>;
    auto coeffs
      = detail::get_main_coeffs (freq, q, db, t_spl, detail::lshelf_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = coeffs.k * (coeffs.A - (T) 1.);
    c[m2] = (coeffs.A * coeffs.A) - (T) 1.;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> t_spl,
    highshelf_tag)
  {
    using T = vec_value_type_t<V>;
    auto coeffs
      = detail::get_main_coeffs (freq, q, db, t_spl, detail::hshelf_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = coeffs.A * coeffs.A;
    c[m1] = coeffs.k * ((T) 1. - coeffs.A) * coeffs.A;
    c[m2] = (T) 1. - (coeffs.A * coeffs.A);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    xspan<const vec_value_type_t<V>> c, // coeffs (1 set)
    xspan<V>                         s, // states (interleaved, SIMD aligned)
    V                                v0)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    V out = calc<V> (
      s[s1],
      s[s2],
      vec_set<V> (c[m0]),
      vec_set<V> (c[m1]),
      vec_set<V> (c[m2]),
      vec_set<V> (c[a1]),
      vec_set<V> (c[a2]),
      vec_set<V> (c[a3]),
      v0);

    return out;
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    xspan<const V> c, // coeffs (interleaved, SIMD aligned)
    xspan<V>       s, // states (interleaved, SIMD aligned)
    V              v0)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    V out
      = calc<V> (s[s1], s[s2], c[m0], c[m1], c[m2], c[a1], c[a2], c[a3], v0);

    return out;
  }
  //----------------------------------------------------------------------------
private:
  template <class T>
  static T calc (
    T& s1_v,
    T& s2_v,
    T  m0_v,
    T  m1_v,
    T  m2_v,
    T  a1_v,
    T  a2_v,
    T  a3_v,
    T  v0)
  {
    auto tick_r = detail::svf_tick (s1_v, s2_v, a1_v, a2_v, a3_v, v0);
    T    out    = m0_v * v0 + m1_v * tick_r.v1 + m2_v * tick_r.v2;
    s1_v        = tick_r.s1;
    s2_v        = tick_r.s2;
    return out;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

}} // namespace artv::andy
