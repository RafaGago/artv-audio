#pragma once

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/blocks/filters/andy_svf.hpp"

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
  static void repair_unsmoothable_coeffs (crange<double>, crange<const double>)
  {}
  //----------------------------------------------------------------------------
  static void tilt (
    crange<double> co,
    double         freq,
    double         q,
    double         gain_db,
    double         sr)
  {
    assert (co.size() >= n_coeffs);
    andy::svf::low_shelf (co, freq, q, gain_db, sr);
    co.shrink_head (andy::svf::n_coeffs);
    andy::svf::high_shelf (co, freq, q, -gain_db, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               v0)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    v0 = andy::svf::tick (co, st, v0);
    co.shrink_head (andy::svf::n_coeffs);
    st.shrink_head (andy::svf::n_states);
    return andy::svf::tick (co, st, v0);
  }
  //----------------------------------------------------------------------------
  static simd_dbl tick (
    crange<double const>          co, // coeffs
    std::array<crange<double>, 2> st, // state
    std::array<double, 2>         v0s)
  {
    assert (st.size() >= 2);
    assert (co.size() >= n_coeffs);
    assert (st[0].size() >= n_states);
    assert (st[1].size() >= n_states);

    auto v0_simd = andy::svf::tick<16> (co, st, v0s);
    co.shrink_head (andy::svf::n_coeffs);
    st[0].shrink_head (andy::svf::n_states);
    st[1].shrink_head (andy::svf::n_states);
    return andy::svf::tick<16> (co, st, {v0_simd[0], v0_simd[1]});
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
