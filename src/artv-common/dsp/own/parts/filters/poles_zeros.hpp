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
  static void reset_coeffs (crange<V> co, vec_complex<V> zero)
  {
    assert (co.size() >= vec_complex<V>::size);

    vec_store (co.data(), zero);
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
  static vec_complex<V> tick (
    crange<const V> co,
    crange<V>       st,
    vec_complex<V>  x)
  {
    assert (co.size() >= vec_complex<V>::vec_size);
    assert (st.size() >= vec_complex<V>::vec_size);

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
  static void reset_coeffs (crange<V> co, vec_complex<V> pole)
  {
    assert (co.size() >= vec_complex<V>::size);

    vec_store (co.data(), pole);
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
  static vec_complex<V> tick (
    crange<const V> co,
    crange<V>       st,
    vec_complex<V>  x)
  {
    assert (co.size() >= vec_complex<V>::vec_size);
    assert (st.size() >= vec_complex<V>::vec_size);

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
  static void reset_coeffs (crange<V> co, V re_zero)
  {
    assert (co.size() >= n_coeffs);
    co[re] = re_zero;
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
  static V tick (crange<const V> co, crange<V> st, V x)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V ret  = x - co[re] * st[z1];
    st[z1] = x;
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
  static void reset_coeffs (crange<V> co, V re_zero)
  {
    assert (co.size() >= n_coeffs);
    co[re] = re_zero;
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
  static V tick (crange<const V> co, crange<V> st, V x)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V ret  = co[re] * st[y1] + x;
    st[y1] = ret;
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
  static void reset_coeffs (crange<V> co, vec_complex<V> pole)
  {
    cpole::reset_coeffs (co, pole);
    co[ratio] = vec_real (pole) / vec_imag (pole);
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
  static V tick (crange<const V> co, crange<V> st, V x)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    vec_complex<V> y = cpole::tick (co, st, vec_complex<V> {x});
    return vec_real (y) + co[ratio] * vec_imag (y);
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
  static void reset_coeffs (crange<V> co, V re_pole, V re_zero)
  {
    rpole::reset_coeffs (co, re_pole);
    co = co.shrink_head (rpole::n_coeffs);
    rzero::reset_coeffs (co, re_zero);
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
  static V tick (crange<const V> co, crange<V> st, V x)
  {
    V out = rpole::tick (co, st, x);
    co    = co.shrink_head (rpole::n_coeffs);
    st    = st.shrink_head (rpole::n_states);
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
  static void reset_coeffs (crange<V> co, vec_complex<V> pole, V re_zero)
  {
    ccpole_pair::reset_coeffs (co, pole);
    co = co.shrink_head (ccpole_pair::n_coeffs);
    rzero::reset_coeffs (co, re_zero);
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
  static V tick (crange<const V> co, crange<V> st, V x)
  {
    V out = ccpole_pair::tick (co, st, x);
    co    = co.shrink_head (ccpole_pair::n_coeffs);
    st    = st.shrink_head (ccpole_pair::n_states);

    for (uint i = 0; i < 2; ++i) {
      out = rzero::tick (co, st, out);
      // same zero location
      st = st.shrink_head (rzero::n_states);
    }
    return out;
  }
  //----------------------------------------------------------------------------
};
} // namespace artv
