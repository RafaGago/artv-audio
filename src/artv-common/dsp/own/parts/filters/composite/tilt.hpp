#pragma once

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"

namespace artv {
//------------------------------------------------------------------------------
// Not strictly a filter, just a linked low and high shelf. Move to another
// place?
class tilt_eq {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs = 2 * andy::svf::n_coeffs;
  static constexpr uint n_states = 2 * andy::svf::n_states;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           q,
    V                           gain_db,
    vec_value_type_t<V>         sr)
  {
    assert (co.size() >= n_coeffs);
    andy::svf::reset_coeffs (co, freq, q, gain_db, sr, lowshelf_tag {});
    co = co.shrink_head (andy::svf::n_coeffs);
    andy::svf::reset_coeffs (co, freq, q, -gain_db, sr, highshelf_tag {});
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
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag              t)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs);
    assert (st.size() >= (n_states * traits.size));

    in = andy::svf::tick (co, st, in, t);
    co = co.shrink_head (andy::svf::n_coeffs);
    st = st.shrink_head (andy::svf::n_states * traits.size);
    return andy::svf::tick (co, st, in, t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // states (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));
    assert (st.size() >= (n_states * traits.size));

    in = andy::svf::tick (co, st, in);
    co = co.shrink_head (andy::svf::n_coeffs * traits.size);
    st = st.shrink_head (andy::svf::n_states * traits.size);
    return andy::svf::tick (co, st, in);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
