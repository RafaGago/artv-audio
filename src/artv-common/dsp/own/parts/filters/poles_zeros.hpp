#pragma once

#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/simd_complex.hpp"

namespace artv {

// complex one-zero filter -----------------------------------------------------
struct czero {
  //----------------------------------------------------------------------------
  enum coeffs { re, im, n_coeffs };
  enum coeffs_int { n_coeffs_int };
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

    auto y = x - zero * z1;
    vec_store (st.data(), x);
    return y;
  }
  //----------------------------------------------------------------------------
};

// complex one-pole filter -----------------------------------------------------
struct cpole {
  //----------------------------------------------------------------------------
  enum coeffs { re, im, n_coeffs };
  enum coeffs_int { n_coeffs_int };
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

    auto y = pole * y1 + x;
    vec_store (st.data(), y);
    return y;
  }
  //----------------------------------------------------------------------------
};
// real one-zero filter --------------------------------------------------------
struct rzero {
  //----------------------------------------------------------------------------
  enum coeffs { re, n_coeffs };
  enum coeffs_int { n_coeffs_int };
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

    V y    = x - co[re] * st[z1];
    st[z1] = x;
    return y;
  }
  //----------------------------------------------------------------------------
};
// real one-pole filter --------------------------------------------------------
struct rpole {
  //----------------------------------------------------------------------------
  enum coeffs { re, n_coeffs };
  enum coeffs_int { n_coeffs_int };
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

    V y    = co[re] * st[y1] + x;
    st[y1] = y;
    return y;
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
  enum coeffs_int { n_coeffs_int };
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
  enum coeffs_int { n_coeffs_int = rzero::n_coeffs_int + rpole::n_coeffs_int };
  enum state { n_states = rzero::n_states + rpole::n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V re_pole, V re_zero)
  {
    rpole::reset_coeffs (co, re_pole);
    co.cut_head (rpole::n_coeffs);
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
    co.cut_head (rpole::n_coeffs);
    st.cut_head (rpole::n_states);
    return rzero::tick (co, st, out);
  }
  //----------------------------------------------------------------------------
};
// complex conjugate poles pair + two real zeros filter ------------------------
struct ccpole_pair_rzero_pair {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = rzero::n_coeffs + ccpole_pair::n_coeffs };
  enum coeffs_int {
    n_coeffs_int = rzero::n_coeffs_int + ccpole_pair::n_coeffs_int
  };
  enum state { n_states = 2 * rzero::n_states + ccpole_pair::n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, vec_complex<V> pole, V re_zero)
  {
    ccpole_pair::reset_coeffs (co, pole);
    co.cut_head (ccpole_pair::n_coeffs);
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
    co.cut_head (ccpole_pair::n_coeffs);
    st.cut_head (ccpole_pair::n_states);

    for (uint i = 0; i < 2; ++i) {
      out = rzero::tick (co, st, out);
      // same zero location
      st.cut_head (rzero::n_states);
    }
    return out;
  }
  //----------------------------------------------------------------------------
};

// complex zero pair, either real or complex conjugate, so it gets real input
// and output.
struct czero_pair {
  //----------------------------------------------------------------------------
  enum coeffs { c1_re, c1_im, c2_re, c2_im, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { z1_c1_re, z1_c2_re, z1_c2_im, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>      co,
    vec_complex<V> zero1,
    vec_complex<V> zero2)
  {
    assert (co.size() >= vec_complex<V>::size);

    vec_store (&co[c1_re], zero1);
    vec_store (&co[c2_re], zero2);
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

    auto zero1 = vec_load<vec_complex<V>> (&co[c1_re]);
    auto zero2 = vec_load<vec_complex<V>> (&co[c2_re]);
    auto z1_2  = vec_load<vec_complex<V>> (&st[z1_c2_re]);

    auto out1    = x - zero1 * st[z1_c1_re];
    st[z1_c1_re] = x;
    auto out2    = out1 - zero2 * z1_2;
    vec_store (&st[z1_c2_re], out1);
    return vec_real (out2);
  }
  //----------------------------------------------------------------------------
};

// complex zero pair, either real or complex conjugate, so it gets real input
// and output.
struct cpole_pair {
  //----------------------------------------------------------------------------
  enum coeffs { c1_re, c1_im, c2_re, c2_im, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1_c1_re, y1_c1_im, y1_c2_re, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>      co,
    vec_complex<V> pole1,
    vec_complex<V> pole2)
  {
    assert (co.size() >= vec_complex<V>::size);

    vec_store (&co[c1_re], pole1);
    vec_store (&co[c2_re], pole2);
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

    auto pole1 = vec_load<vec_complex<V>> (&co[c1_re]);
    auto pole2 = vec_load<vec_complex<V>> (&co[c2_re]);
    auto y     = vec_load<vec_complex<V>> (&st[y1_c1_re]);

    y = pole1 * y + x;
    vec_store (&st[y1_c1_re], y);
    y            = pole2 * st[y1_c2_re] + y;
    V y_re       = vec_real (y);
    st[y1_c2_re] = y_re;
    return y_re;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct cpole_pair_czero_pair {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = cpole_pair::n_coeffs + czero_pair::n_coeffs };
  enum coeffs_int {
    n_coeffs_int = cpole_pair::n_coeffs_int + czero_pair::n_coeffs_int
  };
  enum state { n_states = cpole_pair::n_states + czero_pair::n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>      co,
    vec_complex<V> pole1,
    vec_complex<V> pole2,
    vec_complex<V> zero1,
    vec_complex<V> zero2)
  {
    assert (co.size() >= vec_complex<V>::size);

    czero_pair::reset_coeffs (co, zero1, zero2);
    co.cut_head (czero_pair::n_coeffs);
    cpole_pair::reset_coeffs (co, pole1, pole2);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    czero_pair::reset_states (st);
    st.cut_head (czero_pair::n_states);
    cpole_pair::reset_states (st);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> co, crange<V> st, V x)
  {
    V out = x;
    out   = czero_pair::tick (co, st, out);
    co.cut_head (czero_pair::n_coeffs);
    st.cut_head (czero_pair::n_states);
    out = cpole_pair::tick (co, st, out);
    return out;
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
