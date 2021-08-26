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
    co = co.shrink_head (andy::svf::n_coeffs);
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
    co = co.shrink_head (andy::svf::n_coeffs);
    st = st.shrink_head (andy::svf::n_states);
    return andy::svf::tick (co, st, v0);
  }
  //----------------------------------------------------------------------------
  static double_x2 tick (
    crange<double const>          co, // coeffs
    std::array<crange<double>, 2> st, // state
    double_x2                     v0)
  {
    assert (st.size() >= 2);
    assert (co.size() >= n_coeffs);
    assert (st[0].size() >= n_states);
    assert (st[1].size() >= n_states);

    v0    = andy::svf::tick (co, st, v0);
    co    = co.shrink_head (andy::svf::n_coeffs);
    st[0] = st[0].shrink_head (andy::svf::n_states);
    st[1] = st[1].shrink_head (andy::svf::n_states);
    return andy::svf::tick (co, st, v0);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
