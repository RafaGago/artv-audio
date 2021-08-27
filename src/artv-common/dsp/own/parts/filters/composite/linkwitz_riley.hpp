#pragma once

#include <cassert>
#include <cmath>

#include <algorithm>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

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
  static constexpr uint n_coeffs            = 2 * onepole::n_coeffs;
  static constexpr uint n_states            = 3 * onepole::n_states;
  static constexpr uint n_correction_states = onepole::n_states;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    auto params = onepole::get_sub_coeffs (freq, sr);
    onepole::init (c, params, lowpass_tag {});
    onepole::init (
      c.shrink_head (onepole::n_coeffs * traits.size), params, allpass_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
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
  static lr_crossover_out<V> tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (s.size() >= (n_states * traits.size));

    lr_crossover_out<V> ret;
    // LP
    ret.lp = onepole::tick (c, s, in);
    ret.lp = onepole::tick (
      c, s.shrink_head (onepole::n_states * traits.size), ret.lp);
    // HP = Allpass - LP
    ret.hp = onepole::tick (
      c.shrink_head (onepole::n_coeffs * traits.size),
      s.shrink_head (2 * onepole::n_states * traits.size),
      in);
    ret.hp -= ret.lp;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V apply_correction (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       extern_s, // coeffs "
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (extern_s.size() >= (n_correction_states * traits.size));

    return onepole::tick (
      c.shrink_head (onepole::n_coeffs * traits.size), extern_s, in);
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
  static constexpr uint n_coeffs = svf_lp_ap::n_coeffs;
  static constexpr uint n_states = svf_lp_ap::n_states + svf_lp::n_states;
  static constexpr uint n_correction_states = svf_lp_ap::n_states;

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto        traits = vec_traits<V>();
    static constexpr auto qlist  = butterworth_2p_cascade_q_list::cget<2>();

    assert (c.size() >= (n_coeffs * traits.size));
    svf_lp_ap::init (c, freq, vec_set<V> (qlist[0]), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
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
  static lr_crossover_out<V> tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (s.size() >= (n_states * traits.size));
    lr_crossover_out<V> ret;
    // LP
    auto multimode = svf_lp_ap::tick (c, s, in);
    ret.lp         = svf_lp::tick (
      c, s.shrink_head (svf_lp_ap::n_states * traits.size), multimode[0]);
    // HP = Allpass - LP
    ret.hp = multimode[1] - ret.lp;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V apply_correction (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       extern_s, // coeffs "
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (extern_s.size() >= (n_correction_states * traits.size));

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
  static constexpr uint n_coeffs = 2 * svf_lp_ap::n_coeffs;
  static constexpr uint n_states
    = (2 * svf_lp_ap::n_states) + (3 * svf_lp::n_states);
  static constexpr uint n_correction_states = 2 * svf_lp_ap::n_states;

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto        traits = vec_traits<V>();
    static constexpr auto qlist  = butterworth_2p_cascade_q_list::cget<4>();

    assert (c.size() >= (n_coeffs * traits.size));

    svf_lp_ap::init (c, freq, vec_set<V> (qlist[0]), sr);
    svf_lp_ap::init (
      c.shrink_head (svf_lp_ap::n_coeffs * traits.size),
      freq,
      vec_set<V> (qlist[1]),
      sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
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
  static lr_crossover_out<V> tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 in)
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

    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (s.size() >= (n_states * traits.size));

    lr_crossover_out<V> ret;

    // First Butterworth LP + AP
    auto      coeffs1 = c;
    auto      coeffs2 = c.shrink_head (svf_lp_ap::n_coeffs * traits.size);
    crange<T> s_ptr   = s;

    auto multimode = svf_lp_ap::tick (coeffs1, s_ptr, in);
    s_ptr          = s_ptr.shrink_head (svf_lp_ap::n_states * traits.size);
    V ap           = svf_lp_ap::tick (coeffs2, s_ptr, multimode[1])[1];
    s_ptr          = s_ptr.shrink_head (svf_lp_ap::n_states * traits.size);
    ret.lp         = svf_lp::tick (coeffs2, s_ptr, multimode[0]);
    s_ptr          = s_ptr.shrink_head (svf_lp::n_states * traits.size);

    // Second Butterworth LP
    ret.lp = svf_lp::tick (coeffs1, s_ptr, ret.lp);
    s_ptr  = s_ptr.shrink_head (svf_lp::n_states * traits.size);
    ret.lp = svf_lp::tick (coeffs2, s_ptr, ret.lp);

    // HP = Allpass - LP
    ret.hp = ap - ret.lp;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V apply_correction (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       extern_s, // coeffs "
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (extern_s.size() >= (n_correction_states * traits.size));

    V r = svf_lp_ap::tick (c, extern_s, in)[1];
    return svf_lp_ap::tick (
      c.shrink_head (svf_lp_ap::n_coeffs * traits.size),
      extern_s.shrink_head (svf_lp_ap::n_states * traits.size),
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

  static constexpr uint n_states
    = std::max (lr2::n_states, std::max (lr4::n_states, lr8::n_states));

  static constexpr uint n_correction_states = std::max (
    lr2::n_correction_states,
    std::max (lr4::n_correction_states, lr8::n_correction_states));
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr,
    uint                        order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    switch (order) {
    case 2:
      lr2::init (c, freq, sr);
      break;
    case 4:
      lr4::init (c, freq, sr);
      break;
    case 8:
      lr8::init (c, freq, sr);
      break;
    default:
      assert (false);
      break;
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
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
  static lr_crossover_out<V> tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (s.size() >= (n_states * traits.size));

    switch (order) {
    case 2:
      return lr2::tick (c, s, in);
      break;
    case 4:
      return lr4::tick (c, s, in);
      break;
    case 8:
      return lr8::tick (c, s, in);
      break;
    default:
      assert (false);
      return lr_crossover_out<V> {};
      break;
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V apply_correction (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>> extern_s, // states (interleaved, SIMD aligned)
    V                           in,
    uint                        order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    assert (extern_s.size() >= (n_correction_states * traits.size));

    switch (order) {
    case 2:
      return lr2::apply_correction (c, extern_s, in);
      break;
    case 4:
      return lr4::apply_correction (c, extern_s, in);
      break;
    case 8:
      return lr8::apply_correction (c, extern_s, in);
      break;
    default:
      assert (false);
      return V {};
      break;
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
