#pragma once

#include <cassert>
#include <cmath>

#include <algorithm>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"

namespace artv {

template <class T>
struct lr_crossover_out {
  T lp;
  T hp;
};

template <uint order>
class linkwitz_riley_stage;

//------------------------------------------------------------------------------
// A single crossover split point. A building block for doing crossovers.
//------------------------------------------------------------------------------
template <>
class linkwitz_riley_stage<2> {
public:
  //----------------------------------------------------------------------------
  // uses allpass + lowpass to use the same code everywhere for the corrector,
  // informed-guessing that tighter code is of more value than saving 2
  // subtractions.
  using onepole_type = onepole<lowpass_tag, allpass_tag>;
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs            = 1 * onepole_type::n_coeffs;
  static constexpr uint n_coeffs_int        = 1 * onepole_type::n_coeffs_int;
  static constexpr uint n_states            = 2 * onepole_type::n_states;
  static constexpr uint n_correction_states = onepole_type::n_states;
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> c, V freq, vec_value_type_t<V> t_spl)
  {
    assert (c.size() >= n_coeffs);
    onepole_type::reset_coeffs (c, freq, t_spl);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static lr_crossover_out<V> tick (
    xspan<V const> c, // coeffs (interleaved, SIMD aligned)
    xspan<V>       s, // states (interleaved, SIMD aligned)
    V              in)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    lr_crossover_out<V> ret;
    // LP + AP. lowpasses on index 0, allpasses on 1.
    auto stage1 = onepole_type::tick (c, s, in);
    auto stage2
      = onepole_type::tick (c, s.advanced (onepole_type::n_states), stage1[0]);

    ret.lp = stage2[0];
    ret.hp = stage1[1] - ret.lp; // HP = Allpass - LP
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V apply_correction (
    xspan<V const> c, // coeffs (interleaved, SIMD aligned)
    xspan<V>       extern_s, // coeffs "
    V              in)
  {
    assert (c.size() >= n_coeffs);
    assert (extern_s.size() >= n_correction_states);

    return onepole_type::tick (c, extern_s, in)[1]; // Return allpass out.
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// A single crossover split point. A building block for doing crossovers.
//------------------------------------------------------------------------------
template <>
class linkwitz_riley_stage<4> {
private:
  //----------------------------------------------------------------------------
  using svf_lp_ap = andy::svf_multimode<lowpass_tag, allpass_tag>;
  using svf_lp    = andy::svf_lowpass;

public:
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = svf_lp_ap::n_coeffs;
  static constexpr uint n_coeffs_int = svf_lp_ap::n_coeffs_int;
  static constexpr uint n_states     = svf_lp_ap::n_states + svf_lp::n_states;
  static constexpr uint n_correction_states = svf_lp_ap::n_states;

  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> c, V freq, vec_value_type_t<V> t_spl)
  {
    static constexpr auto qlist = butterworth_2p_cascade_q_list::cget<2>();

    assert (c.size() >= n_coeffs);
    svf_lp_ap::reset_coeffs (c, freq, vec_set<V> (qlist[0]), t_spl);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static lr_crossover_out<V> tick (
    xspan<V const> c, // coeffs (interleaved, SIMD aligned)
    xspan<V>       s, // states (interleaved, SIMD aligned)
    V              in)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);
    lr_crossover_out<V> ret;
    // LP
    auto multimode = svf_lp_ap::tick (c, s, in);
    ret.lp = svf_lp::tick (c, s.advanced (svf_lp_ap::n_states), multimode[0]);
    // HP = Allpass - LP
    ret.hp = multimode[1] - ret.lp;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V apply_correction (
    xspan<V const> c, // coeffs (interleaved, SIMD aligned)
    xspan<V>       extern_s, // coeffs "
    V              in)
  {
    assert (c.size() >= n_coeffs);
    assert (extern_s.size() >= n_correction_states);

    return svf_lp_ap::tick (c, extern_s, in)[1];
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// A single crossover split point. A building block for doing crossovers.
//------------------------------------------------------------------------------
template <>
class linkwitz_riley_stage<8> {
private:
  //----------------------------------------------------------------------------
  using svf_lp_ap = andy::svf_multimode<lowpass_tag, allpass_tag>;
  using svf_lp    = andy::svf_lowpass;

public:
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 2 * svf_lp_ap::n_coeffs;
  static constexpr uint n_coeffs_int = 2 * svf_lp_ap::n_coeffs_int;
  static constexpr uint n_states
    = (2 * svf_lp_ap::n_states) + (3 * svf_lp::n_states);
  static constexpr uint n_correction_states = 2 * svf_lp_ap::n_states;

  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> c, V freq, vec_value_type_t<V> t_spl)
  {
    static constexpr auto qlist = butterworth_2p_cascade_q_list::cget<4>();

    assert (c.size() >= n_coeffs);

    svf_lp_ap::reset_coeffs (c, freq, vec_set<V> (qlist[0]), t_spl);
    svf_lp_ap::reset_coeffs (
      c.advanced (svf_lp_ap::n_coeffs), freq, vec_set<V> (qlist[1]), t_spl);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static lr_crossover_out<V> tick (
    xspan<V const> c, // coeffs (interleaved, SIMD aligned)
    xspan<V>       s, // states (interleaved, SIMD aligned)
    V              in)
  {
    // Look at:
    // https://www.kvraudio.com/forum/viewtopic.php?t=557465
    //
    // Notice that R1 and R2 is the inverse of the Q's of each butterworth
    // filter. so "R1=1/0.54119610" and "R2=1/1.3065630".
    //
    // Copying from the thread the signal chain is:
    //
    // ---LR8---
    // SVF2pole1(in,R1 -> hp2,bp2,lp2)
    // SVF2pole2(lp2,R2 -> lp2hp2,lp2bp2,lp4)
    // SVF2pole3(lp4,R1 -> lp6)
    // SVF2pole4(lp6,R2 -> lp8)
    // ap4 = m0 hp2 + m1 bp2 + m3 lp2hp2 + m4 lp2bp2 + m5 lp4
    // hp8 = ap4 - lp8
    //
    // Where:
    //
    // {m0 = 1, m1 = -2 (R1 + 2 R2), m2 = 0, m3 = 1 + 8 R2 (R1 + R2), m4 = 2 R2,
    //  m5 = 1}
    //
    // The number of total operations stays more or less the same though, and
    // the corrective allpass can run in parallel once the first LP is
    // processed, so if going for this optimization, measure performance.

    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    lr_crossover_out<V> ret;

    // First Butterworth LP + AP
    auto coeffs1 = c;
    auto coeffs2 = c.advanced (svf_lp_ap::n_coeffs);
    auto s_ptr   = s;

    auto multimode = svf_lp_ap::tick (coeffs1, s_ptr, in);
    s_ptr.cut_head (svf_lp_ap::n_states);

    V ap = svf_lp_ap::tick (coeffs2, s_ptr, multimode[1])[1];
    s_ptr.cut_head (svf_lp_ap::n_states);

    ret.lp = svf_lp::tick (coeffs2, s_ptr, multimode[0]);
    s_ptr.cut_head (svf_lp::n_states);

    // Second Butterworth LP
    ret.lp = svf_lp::tick (coeffs1, s_ptr, ret.lp);
    s_ptr.cut_head (svf_lp::n_states);

    ret.lp = svf_lp::tick (coeffs2, s_ptr, ret.lp);

    // HP = Allpass - LP
    ret.hp = ap - ret.lp;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V apply_correction (
    xspan<V const> c, // coeffs (interleaved, SIMD aligned)
    xspan<V>       extern_s, // coeffs "
    V              in)
  {
    assert (c.size() >= n_coeffs);
    assert (extern_s.size() >= n_correction_states);

    V r = svf_lp_ap::tick (c, extern_s, in)[1];
    return svf_lp_ap::tick (
      c.advanced (svf_lp_ap::n_coeffs),
      extern_s.advanced (svf_lp_ap::n_states),
      r)[1];
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// A single crossover split point. A building block for doing crossovers.
//------------------------------------------------------------------------------
class linkwitz_riley_stage_any_order {
public:
  //----------------------------------------------------------------------------
  using lr2 = linkwitz_riley_stage<2>;
  using lr4 = linkwitz_riley_stage<4>;
  using lr8 = linkwitz_riley_stage<8>;
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs
    = std::max (lr2::n_coeffs, std::max (lr4::n_coeffs, lr8::n_coeffs));

  static constexpr uint n_coeffs_int = std::max (
    lr2::n_coeffs_int,
    std::max (lr4::n_coeffs_int, lr8::n_coeffs_int));

  static constexpr uint n_states
    = std::max (lr2::n_states, std::max (lr4::n_states, lr8::n_states));

  static constexpr uint n_correction_states = std::max (
    lr2::n_correction_states,
    std::max (lr4::n_correction_states, lr8::n_correction_states));
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    vec_value_type_t<V> t_spl,
    uint                order)
  {
    assert (c.size() >= n_coeffs);

    switch (order) {
    case 2:
      lr2::reset_coeffs (c, freq, t_spl);
      break;
    case 4:
      lr4::reset_coeffs (c, freq, t_spl);
      break;
    case 8:
      lr8::reset_coeffs (c, freq, t_spl);
      break;
    default:
      assert (false);
      break;
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<vec_value_type_t<V>> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static lr_crossover_out<V> tick (
    xspan<V const> c,
    xspan<V>       s,
    V              in,
    uint           order)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    switch (order) {
    case 2:
      return lr2::tick (c, s, in);
    case 4:
      return lr4::tick (c, s, in);
    case 8:
      return lr8::tick (c, s, in);
    default:
      assert (false);
      return lr_crossover_out<V> {};
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V apply_correction (
    xspan<V const> c,
    xspan<V>       extern_s,
    V              in,
    uint           order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= n_coeffs);
    assert (extern_s.size() >= n_correction_states);

    switch (order) {
    case 2:
      return lr2::apply_correction (c, extern_s, in);
    case 4:
      return lr4::apply_correction (c, extern_s, in);
    case 8:
      return lr8::apply_correction (c, extern_s, in);
    default:
      assert (false);
      return V {};
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
