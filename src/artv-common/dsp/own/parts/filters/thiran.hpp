#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/traits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

//------------------------------------------------------------------------------
namespace artv {

template <uint order>
class thiran;

// direct form 1 for more robust coefficient changes
template <>
class thiran<1> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { a, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1, x1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co,
    V                   freq,
    vec_value_type_t<V> t_spl,
    quality_tag<0>) // no frequency prewarp
  {
    // as allpass filter (negative coeff)
    using T = vec_value_type_t<V>;
    auto d  = (T) 2 * freq * t_spl;
    co[a]   = (d - (T) 1) / (d + (T) 1);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V freq, vec_value_type_t<V> t_spl)
  {
    // as allpass filter (negative coeff)
    using T = vec_value_type_t<V>;
    auto d  = vec_tan ((T) M_PI * freq * t_spl);
    co[a]   = (d - (T) 1) / (d + (T) 1);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V fractional)
  {
    // as interpolator
    using T = vec_value_type_t<V>;
    // D parameter between 0.418 and 1.418
    V d = fractional + (T) 0.418;
#if 1
    co[a] = ((T) 1 - d) / ((T) 1 + d);
#else
    // See https://dafx09.como.polimi.it/proceedings/papers/paper_72.pdf
    // chapter 6.
    V v = ((T) 1. - d) * (T) (1. / 2.);

    V v_p2 = v * v;
    V v_p4 = v_p2 * v_p2;
    V v_p8 = v_p4 * v_p4;

    co[a] = v;
    co[a] *= ((T) 1 + v);
    co[a] *= ((T) 1 + v_p2);
    co[a] *= ((T) 1 + v_p4);
    co[a] *= ((T) 1 + v_p8);
#endif
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
    xspan<V const> co, // coeffs (interleaved, SIMD aligned)
    xspan<V>       st, // states (interleaved, SIMD aligned)
    V              in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    V y    = in * co[a] + st[x1] - co[a] * st[y1];
    st[y1] = y;
    st[x1] = in;
    return y;
  }
  //----------------------------------------------------------------------------
};
// direct form 1 for more robust coefficient changes
template <>
class thiran<2> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1, y2, x1, x2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V fractional)
  {
    // as interpolator
    using T = vec_value_type_t<V>;
    V d     = fractional + (T) 1.403;
#if 1
    co[a1] = -(d - (T) 2) / (d + (T) 1);
    co[a2] = ((d - (T) 1) * (d - (T) 2)) / ((d + (T) 1) * (d + (T) 2));
#else
    // See https://dafx09.como.polimi.it/proceedings/papers/paper_72.pdf
    // chapter 6.
    // These measured worse than the divisions on a Ryzen7-5800x running 8x16
    // delay lines.
    V v1 = ((T) 3. - d) * (T) (1. / 4.);
    V v2 = ((T) 2. - d) * (T) (1. / 4.);

    V v1_p2 = v1 * v1;
    V v1_p4 = v1_p2 * v1_p2;
    V v1_p8 = v1_p4 * v1_p4;

    co[a1] = (T) -2 * ((T) 0.25 - v1);
    co[a1] *= ((T) 1 + v1);
    co[a1] *= ((T) 1 + v1_p2);
    co[a1] *= ((T) 1 + v1_p4);
    co[a1] *= ((T) 1 + v1_p8);

    V v2_p2 = v2 * v2;
    V v2_p4 = v2_p2 * v2_p2;
    V v2_p8 = v2_p4 * v2_p4;

    co[a2] = co[a1] * (T) -0.5 * ((T) 0.25 - v2);
    co[a2] *= ((T) 1 + v2);
    co[a2] *= ((T) 1 + v2_p2);
    co[a2] *= ((T) 1 + v2_p4);
    co[a2] *= ((T) 1 + v2_p8);
#endif
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
    xspan<V const> co, // coeffs (interleaved, SIMD aligned)
    xspan<V>       st, // states (interleaved, SIMD aligned)
    V              in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    V y = co[a2] * in + co[a1] * st[x1] + st[x2];
    y -= co[a1] * st[y1] + co[a2] * st[y2];
    st[x2] = st[x1];
    st[x1] = in;
    st[y2] = st[y1];
    st[y1] = y;
    return y;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

} // namespace artv
