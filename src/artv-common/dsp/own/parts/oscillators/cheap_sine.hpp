#pragma once

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// https://www.musicdsp.org/en/latest/Synthesis/10-fast-sine-and-cosine-calculation.html
//------------------------------------------------------------------------------
class cheap_sin_osc {
public:
  //----------------------------------------------------------------------------
  using value_type = float;
  //----------------------------------------------------------------------------
  enum coeffs { a, n_coeffs };
  enum state { z0, z1, n_states };
  //----------------------------------------------------------------------------
  static void calc_coefs (crange<value_type> c, float freq, float srate)
  {
    assert (c.size() >= n_coeffs);
    c[a] = 2.f * (float) sin (M_PI * freq / srate);
  }
  //----------------------------------------------------------------------------
  static void reset_states (crange<value_type> st) { reset_states (st, 1.0f); }
  //----------------------------------------------------------------------------
  static void reset_states (crange<value_type> st, float bipolar_ampl)
  {
    assert (st.size() >= n_states);

    st[z0] = bipolar_ampl; // -1 to +1 oscillation as default
    st[z1] = 0.f;
  }
  //----------------------------------------------------------------------------
  static float tick (crange<const value_type> co, crange<float> st)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    st[z0] = st[z0] - co[a] * st[z1];
    st[z1] = st[z1] + co[a] * st[z0];
    return st[z0];
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
