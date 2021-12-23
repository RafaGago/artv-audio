#pragma once

#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/simd_complex.hpp"

namespace artv {

// complex one-zero filter -----------------------------------------------------
struct czero {
  //----------------------------------------------------------------------------
  enum coeffs { re, im, n_coeffs };
  enum state { z1_re, z1_im, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, vec_complex<V> zero)
  {
    assert (co.size() >= vec_complex<V>::size);

    vec_store (co.data(), zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static vec_complex<V> tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    vec_complex<V>                    x)
  {
    assert (co.size() >= vec_complex<V>::size);
    assert (st.size() >= vec_complex<V>::size);

    auto zero = vec_load<vec_complex<V>> (co.data());
    auto z1   = vec_load<vec_complex<V>> (st.data());

    auto ret = x - zero * z1;
    vec_store (st.data(), x);
    return ret;
  }
  //----------------------------------------------------------------------------
};
// complex one-pole filter -----------------------------------------------------
struct cpole {
  //----------------------------------------------------------------------------
  enum coeffs { re, im, n_coeffs };
  enum state { y1_re, y1_im, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, vec_complex<V> pole)
  {
    assert (co.size() >= vec_complex<V>::size);

    vec_store (co.data(), pole);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static vec_complex<V> tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    vec_complex<V>                    x)
  {
    assert (co.size() >= vec_complex<V>::size);
    assert (st.size() >= vec_complex<V>::size);

    auto pole = vec_load<vec_complex<V>> (co.data());
    auto y1   = vec_load<vec_complex<V>> (st.data());

    auto ret = pole * y1 + x;
    vec_store (st.data(), ret);
    return ret;
  }
  //----------------------------------------------------------------------------
};
// real one-zero filter --------------------------------------------------------
struct rzero {
  //----------------------------------------------------------------------------
  enum coeffs { re, n_coeffs };
  enum state { z1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, V re_zero)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    vec_store (&co[re * traits.size], re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    V re_v = vec_load<V> (&co[re * traits.size]);
    V z1_v = vec_load<V> (&st[z1 * traits.size]);

    V ret = x - re_v * z1_v;
    vec_store (&st[z1 * traits.size], x);
    return ret;
  }
  //----------------------------------------------------------------------------
};
// real one-pole filter --------------------------------------------------------
struct rpole {
  //----------------------------------------------------------------------------
  enum coeffs { re, n_coeffs };
  enum state { y1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, V re_zero)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    vec_store (&co[re * traits.size], re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    V re_v = vec_load<V> (&co[re * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);

    V ret = re_v * y1_v + x;
    vec_store (&st[y1 * traits.size], ret);
    return ret;
  }
  //----------------------------------------------------------------------------
};
// complex conjugate pole pair filter ------------------------------------------
// Based on paper from Martin Vicanek "A New Reverse IIR Algorithm".
// https://www.vicanek.de/articles.htm
// file:///home/s0001192/Downloads/ReverseIIR.pdf
struct ccpole_pair {
  //----------------------------------------------------------------------------
  enum coeffs { ratio = cpole::n_coeffs, n_coeffs };
  enum state { n_states = cpole::n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, vec_complex<V> pole)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    cpole::reset_coeffs (co, pole);
    vec_store (&co[ratio * traits.size], vec_real (pole) / vec_imag (pole));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * n_states));

    vec_complex<V> y = cpole::tick (co, st, vec_complex<V> {x});
    return vec_real (y) + vec_load<V> (&co[ratio * traits.size]) * vec_imag (y);
  }
  //----------------------------------------------------------------------------
};
// one real pole one real zero filter ------------------------------------------
struct rpole_rzero {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = rzero::n_coeffs + rpole::n_coeffs };
  enum state { n_states = rzero::n_states + rpole::n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           re_pole,
    V                           re_zero)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    rpole::reset_coeffs (co, re_pole);
    co = co.shrink_head (rpole::n_coeffs * traits.size);
    rzero::reset_coeffs (co, re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    V out = rpole::tick (co, st, x);
    co    = co.shrink_head (rpole::n_coeffs * traits.size);
    st    = st.shrink_head (rpole::n_states * traits.size);
    return rzero::tick (co, st, out);
  }
  //----------------------------------------------------------------------------
};
// complex conjugate poles pair + two real zeros filter ------------------------
struct ccpole_pair_rzero_pair {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = rzero::n_coeffs + ccpole_pair::n_coeffs };
  enum state { n_states = 2 * rzero::n_states + ccpole_pair::n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    vec_complex<V>              pole,
    V                           re_zero)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    ccpole_pair::reset_coeffs (co, pole);
    co = co.shrink_head (ccpole_pair::n_coeffs * traits.size);
    rzero::reset_coeffs (co, re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 x)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    V out = ccpole_pair::tick (co, st, x);
    co    = co.shrink_head (ccpole_pair::n_coeffs * traits.size);
    st    = st.shrink_head (ccpole_pair::n_states * traits.size);

    for (uint i = 0; i < 2; ++i) {
      out = rzero::tick (co, st, out);
      // same zero location
      st = st.shrink_head (rzero::n_states * traits.size);
    }
    return out;
  }
  //----------------------------------------------------------------------------
};
} // namespace artv
