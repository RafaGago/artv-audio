#pragma once

#include <array>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/filters/zdf.hpp"
#include "artv-common/dsp/own/parts/traits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

//------------------------------------------------------------------------------
namespace artv {

struct onepole_smoother {
  //----------------------------------------------------------------------------
  enum coeffs { b1, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c, V freq, vec_value_type_t<V> t_spl)
  {
    using T              = vec_value_type_t<V>;
    constexpr auto pi_x2 = (T) (2. * M_PI);

    c[b1] = vec_exp (-pi_x2 * freq * t_spl);
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
  static V tick (
    crange<const vec_value_type_t<V>> c, // coeffs (1 set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T = vec_value_type_t<V>;

    assert (st.size() >= n_states);
    assert (c.size() >= n_coeffs);

    V a0_v = vec_set<V> (((T) 1.) - c[b1]);
    V b1_v = vec_set<V> (c[b1]);

    st[y1] = (in * a0_v) + (st[y1] * b1_v);
    return st[y1];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> c, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               in)
  {
    assert (st.size() >= n_states);
    assert (c.size() >= n_coeffs);

    st[y1] = (in * (1. - c[b1])) + (st[y1] * c[b1]);
    return st[y1];
  }
  //----------------------------------------------------------------------------
};
namespace detail {
//------------------------------------------------------------------------------
// Chapter 3.10 THE ART OF VA FILTER DESIGN
//
// Vadim Zavalishin
// https://www.native-instruments.com/fileadmin/ni_media/downloads/pdf/VAFilterDesign_2.1.0.pdf
//
// Notice that this multimode only supports naive shelves, otherwise the cutoff
// frequency for the other outputs would have to be corrected, invalidating the
// multimode.
//------------------------------------------------------------------------------
template <class Derived, class... Tags>
struct onepole {
  //----------------------------------------------------------------------------
  using enabled_modes = mp11::mp_unique<mp11::mp_list<Tags...>>;

  template <class tag>
  static constexpr bool contains = mp11::mp_contains<enabled_modes, tag>::value;

  static constexpr bool returns_array = mp11::mp_size<enabled_modes>::value > 1;
  static constexpr bool has_shelves
    = contains<lowshelf_naive_tag> || contains<highshelf_naive_tag>;

  static_assert (
    mp11::mp_size<enabled_modes>::value >= 1,
    "One mode needs to be enabled");
  //----------------------------------------------------------------------------
  enum coeffs { G, k, n_coeffs = k + (int) has_shelves };
  enum coeffs_int { n_coeffs_int };
  enum state { s, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  // raw overload
  static void reset_coeffs (crange<V> co, V g)
  {
    static_assert (!has_shelves);
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);
    co[G] = g / ((T) 1.0 + g);
    if constexpr (!std::is_same_v<Derived, void>) {
      Derived::after_reset_coeffs (co, g);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V freq, vec_value_type_t<V> t_spl)
  {
    static_assert (!has_shelves);
    using T = vec_value_type_t<V>;
    reset_coeffs (co, vec_tan ((T) M_PI * freq * t_spl));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    vec_value_type_t<V> t_spl,
    no_prewarp)
  {
    static_assert (!has_shelves);
    using T = vec_value_type_t<V>;
    reset_coeffs (co, (T) M_PI * freq * t_spl);
  }
  //----------------------------------------------------------------------------
  // raw overload
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V g, V k_)
  {
    static_assert (has_shelves);
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);
    co[G] = g / ((T) 1.0 + g);
    co[k] = k_;
    if constexpr (!std::is_same_v<Derived, void>) {
      Derived::after_reset_coeffs (co, g);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   db,
    vec_value_type_t<V> t_spl,
    no_prewarp)
  {
    static_assert (has_shelves);
    using T = vec_value_type_t<V>;
    reset_coeffs (
      co, (T) M_PI * freq * t_spl, vec_exp (db * (T) (M_LN10 / 20.)) - (T) 1);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   db,
    vec_value_type_t<V> t_spl)
  {
    static_assert (has_shelves);
    using T = vec_value_type_t<V>;
    reset_coeffs (
      co,
      vec_tan ((T) M_PI * freq * t_spl),
      vec_exp (db * (T) (M_LN10 / 20.)) - (T) 1);
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
  static auto tick (crange<const V> co, crange<V> st, V in)
  {
    return tick_impl<V, V> (co, st, in);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (crange<const vec_value_type_t<V>> co, crange<V> st, V in)
  {
    return tick_impl<V, vec_value_type_t<V>> (co, st, in);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class V, class VT, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick_impl (
    crange<const VT> co, // coeffs (V builtin type (single set) or V (SIMD))
    crange<V>        st, // states (interleaved, SIMD aligned)
    V                x)
  {
    using T = vec_value_type_t<V>;

    V G_, k_, lp, hp;
    if constexpr (std::is_same_v<VT, V>) {
      G_ = co[G];
      if constexpr (has_shelves) {
        k_ = co[k];
      }
    }
    else if constexpr (std::is_same_v<VT, T>) {
      G_ = vec_set<V> (co[G]);
      if constexpr (has_shelves) {
        k_ = vec_set<V> (co[k]);
      }
    }
    else {
      static_assert (
        sizeof (V) == 0, "T type must be either V or V's builtin type");
    }

    V v   = (x - st[s]) * G_;
    lp    = v + st[s];
    st[s] = lp + v;

    std::array<V, mp11::mp_size<enabled_modes>::value> ret;

    mp_foreach_idx (enabled_modes {}, [&] (auto index, auto mode) {
      using tag = decltype (mode);

      if constexpr (std::is_same_v<tag, lowpass_tag>) {
        ret[index.value] = lp;
      }
      else if constexpr (std::is_same_v<tag, highpass_tag>) {
        // The optimizer should trivially detect hp being calculated multiple
        // times. Avoiding throwing templates to the problem.
        hp               = x - lp;
        ret[index.value] = hp;
      }
      else if constexpr (std::is_same_v<tag, allpass_tag>) {
        // The optimizer should trivially detect hp being calculated multiple
        // times. Avoiding throwing templates to the problem.
        hp               = x - lp;
        ret[index.value] = lp - hp;
      }
      else if constexpr (std::is_same_v<tag, lowshelf_naive_tag>) {
        ret[index.value] = x + k_ * lp;
      }
      else if constexpr (std::is_same_v<tag, highshelf_naive_tag>) {
        // The optimizer should trivially detect hp being calculated multiple
        // times. Avoiding throwing templates to the problem.
        hp               = x - lp;
        ret[index.value] = x + k_ * hp;
      }
    });

    if constexpr (returns_array) {
      return ret;
    }
    else {
      return ret[0];
    }
  }
};

} // namespace detail

template <class... Tags>
using onepole = detail::onepole<void, Tags...>;
//------------------------------------------------------------------------------
using onepole_lowpass  = onepole<lowpass_tag>;
using onepole_highpass = onepole<highpass_tag>;
using onepole_allpass  = onepole<allpass_tag>;
//------------------------------------------------------------------------------
// 1 pole adding the ability to get the S and G parameters for zdf feedback
// calculations. Using a separate class because it needs some extra coefficient
// memory to avoid divisions
template <class... Tags>
class onepole_zdf : private detail::onepole<onepole_zdf<Tags...>, Tags...> {
public:
  using base = detail::onepole<onepole_zdf<Tags...>, Tags...>;

  template <class Tag>
  static constexpr bool contains = base::template contains<Tag>;

  static constexpr bool needs_ap_coeff = contains<allpass_tag>;
  // clang-format off
  static constexpr bool needs_g0_coeff =
    contains<highpass_tag> ||
    contains<lowshelf_naive_tag> ||
    contains<highshelf_naive_tag>;
  // clang-format on

  enum coeffs {
    g0       = base::n_coeffs,
    g_ap     = g0 + (uint) needs_g0_coeff,
    n_coeffs = g_ap + (uint) needs_ap_coeff
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
  static auto tick (crange<const V> co, crange<V> st, zdf::gs_coeffs_tag)
  {
    std::array<V, n_zdf_coeffs> ret;
    tick_impl<V, V> (co, st, ret, zdf::gs_coeffs_tag {});
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> co,
    crange<V>                         st,
    zdf::gs_coeffs_tag)
  {
    std::array<V, n_zdf_coeffs> ret;
    tick_impl<V, vec_value_type_t<V>> (co, st, ret, zdf::gs_coeffs_tag {});
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    crange<const V> co,
    crange<V>       st,
    crange<V>       G_S,
    zdf::gs_coeffs_tag)
  {
    tick_impl<V, V> (co, st, G_S, zdf::gs_coeffs_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    crange<const vec_value_type_t<V>> co,
    crange<V>                         st,
    crange<V>                         G_S,
    zdf::gs_coeffs_tag)
  {
    tick_impl<V, vec_value_type_t<V>> (co, st, G_S, zdf::gs_coeffs_tag {});
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  // return G and S for each mode
  template <class V, class VT, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick_impl (
    crange<const VT> co, // coeffs (V builtin type (single set) or V (SIMD))
    crange<V>        st, // states (interleaved, SIMD aligned)
    crange<V>        G_S,
    zdf::gs_coeffs_tag)
  {
    using T = vec_value_type_t<V>;
    assert (G_S.size() >= n_zdf_coeffs);

    V G_, g0_, k_, g_ap_;
    if constexpr (std::is_same_v<VT, V>) {
      G_ = co[base::G];
      if constexpr (needs_g0_coeff) {
        g0_ = co[g0];
      }
      if constexpr (base::has_shelves) {
        k_ = co[base::k];
      }
      if constexpr (needs_ap_coeff) {
        g_ap_ = co[g_ap];
      }
    }
    else if constexpr (std::is_same_v<VT, T>) {
      G_ = vec_set<V> (co[base::G]);
      if constexpr (needs_g0_coeff) {
        g0_ = vec_set<V> (co[g0]);
      }
      if constexpr (base::has_shelves) {
        k_ = vec_set<V> (co[base::k]);
      }
      if constexpr (needs_ap_coeff) {
        g_ap_ = vec_set<V> (co[g_ap]);
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

        if constexpr (std::is_same_v<tag, lowpass_tag>) {
          //          g⋅(s₁ + x)
          // LP_1_y = ───────────
          //             g + 1
          //            g
          // LP_1_G = ─────
          //          g + 1
          //          g⋅s₁
          // LP_1_S = ──────
          //          g + 1
          G_S[i + 0] = G_; // G
          G_S[i + 1] = st[base::s] * G_; // S
        }
        else if constexpr (std::is_same_v<tag, highpass_tag>) {
          //          g⋅s₁ + x
          // HP_1_y = ────────
          //           g + 1
          //            1
          // HP_1_G = ─────
          //          g + 1
          //           -g⋅s₁
          // HP_1_S = ─────
          //          g + 1
          G_S[i + 0] = g0_; // G
          G_S[i + 1] = -st[base::s] * G_; // S
        }
        else if constexpr (std::is_same_v<tag, allpass_tag>) {
          //          -2⋅g⋅s₁ + g⋅x - x
          // AP_1_y = ─────────────────
          //                g + 1
          //          g - 1
          // AP_1_G = ─────
          //          g + 1
          //          -2⋅g⋅s₁
          // AP_1_S = ────────
          //           g + 1
          G_S[i + 0] = g_ap_;
          G_S[i + 1] = (T) 2. * st[base::s] * G_; // S
        }
        else if constexpr (std::is_same_v<tag, lowshelf_naive_tag>) {
          //          -g⋅k⋅(s₁ - x) + x⋅(g + 1)
          // LS_1_y = ─────────────────────────
          //                    g + 1
          //          g⋅k + g + 1
          // LS_1_G = ───────────
          //             g + 1
          //          g⋅k⋅s₁
          // LS_1_S = ────────
          //           g + 1
          G_S[i + 0] = G_ * ((T) 1 + k_) + g0_; // G
          G_S[i + 1] = st[base::s] * k_ * G_; // S
        }
        else if constexpr (std::is_same_v<tag, highshelf_naive_tag>) {
          // y = x + k * yHP
          //          g⋅k⋅s₁ + g⋅x + k⋅x + x
          // HS_1_y = ──────────────────────
          //                  g + 1
          //          g + k + 1
          // HS_1_G = ─────────
          //            g + 1
          //          -g⋅k⋅s₁
          // HS_1_S = ──────
          //          g + 1
          G_S[i + 0] = g0_ * ((T) 1 + k_) + G_; // G
          G_S[i + 1] = -st[base::s] * k_ * G_; // S
        }
      });
  }
  //----------------------------------------------------------------------------
  friend class detail::onepole<onepole_zdf<Tags...>, Tags...>;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static inline void after_reset_coeffs (crange<V> co, V g)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);
    // this div is trivial to remove for the optizer, as g/(1+g) happens just
    // before.
    if constexpr (needs_g0_coeff) {
      co[g0] = (T) 1.0 / ((T) 1.0 + g);
    }
    if constexpr (needs_ap_coeff) {
      co[g_ap] = (g - (T) 1.0) * co[g0];
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// When possible prfer the TPT variant, it only has 1 coeff and 1 state.
struct onepole_tdf2 {
  //----------------------------------------------------------------------------
  enum coeffs { b0, b1, a1, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { s1, n_states };

  template <class T>
  struct sub_coeffs {
    T w;
    T n;
  };
  //----------------------------------------------------------------------------
  template <class T, class U>
  static sub_coeffs<T> get_sub_coeffs (T freq, U t_spl)
  {
    sub_coeffs<T> r;
    if constexpr (is_vec_v<T>) {
      r.w = vec_tan ((U) M_PI * freq * t_spl);
    }
    else {
      r.w = tan (M_PI * freq * t_spl);
    }
    r.n = (U) 1. / ((U) 1. + r.w);
    return r;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, sub_coeffs<V> wn, lowpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    co[b0] = wn.w * wn.n;
    co[b1] = co[b0];
    co[a1] = wn.n * (wn.w - (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, sub_coeffs<V> wn, highpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    co[b0] = wn.n;
    co[b1] = -wn.n;
    co[a1] = wn.n * (wn.w - (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, sub_coeffs<V> wn, allpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    // Found empirically
    co[b0] = wn.n * (wn.w - (T) 1.);
    co[b1] = vec_set<V> ((T) 1.);
    co[a1] = co[b0];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    vec_value_type_t<V> t_spl,
    lowpass_tag         t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, t_spl), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    vec_value_type_t<V> t_spl,
    highpass_tag        t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, t_spl), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    vec_value_type_t<V> t_spl,
    allpass_tag         t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, t_spl), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V> co,
    V         d, // 0 to 1
    thiran_tag)
  {
    using T = vec_value_type_t<V>;

    // http://users.spa.aalto.fi/vpv/publications/vesan_vaitos/ch3_pt3_allpass.pdf
    d += (T) 0.418;
    co[a1] = ((T) 1 - d) / ((T) 1 + d);
    co[b0] = co[a1];
    co[b1] = vec_set<V> ((T) 0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V a1v, V b0v, V b1v, raw_tag)
  {
    co[a1] = a1v;
    co[b0] = b0v;
    co[b1] = b1v;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V b0v = vec_set<V> (co[b0]);
    V b1v = vec_set<V> (co[b1]);
    V a1v = vec_set<V> (co[a1]);

    auto out = in * b0v + st[s1];
    st[s1]   = in * b1v - out * a1v;
    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    auto out = in * co[b0] + st[s1];
    st[s1]   = in * co[b1] - out * co[a1];
    return out;
  }
};
//------------------------------------------------------------------------------
} // namespace artv
