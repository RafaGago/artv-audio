// Ported from https://github.com/Rcomian/jsfx-vcv
// commit sha: 65bfa7896480f3de7d424e42a6ef2b3737b43a20
#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
// slider1:0<-12, 12, 0.1>Input
// slider2:0<0, 10, 0.1>Low Contour
// slider3:0<0, 10, 0.1>Process
// slider4:2<0, 10, 0.1>CV
// slider5:0<-12, 12, 0.1>Output
// slider6:1<0, 1, {Off,On})>Noise Floor

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/third_party/jsfx_engine/jsfx_engine.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv { namespace liteon {

struct bbe {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::exciter;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
private:
#if 0
  /* "mempool" emulates the JSFX per instance memory pool.
  It can be removed if unused. On most cases the length of this array can be
  reduced. It is set to a very conservative (high) value matching what JSFX
  provides. */
  double mempool[8 * 1024 * 1024] = {};
#endif
  //----------------------------------------------------------------------------
  // definitions for environment function calls
  static double jsfx_rand (double maxv = 1.)
  {
    return jsfx_engine::rand (maxv);
  }
  //----------------------------------------------------------------------------
  // stubs for sliders
  //----------------------------------------------------------------------------
public:
#if 0
  void slider()
  {
    g_in        = std::pow (2., get_slider_slider1() * r_6);
    g_lp        = get_slider_slider2() * 0.5;
    g_hp        = get_slider_slider3() * 0.5;
    cv_k        = get_slider_slider4();
    g_out       = std::pow (2., get_slider_slider5() * r_6);
    noise_floor = get_slider_slider6();
  }
#endif
  //----------------------------------------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float v)
  {
    g_in  = std::pow (2., v * r_6);
    g_out = std::pow (2., -v * r_6);
  }

  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -24.f, 24.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct low_contour_tag {};
  void set (low_contour_tag, float v)
  {
    g_lp = v / 20.; // = 0 to 5 range
  }

  static constexpr auto get_parameter (low_contour_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct process_tag {};
  void set (process_tag, float v)
  {
    g_hp = v / 20.; // = 0 to 5 range
  }

  static constexpr auto get_parameter (process_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct cv_tag {};
  void set (cv_tag, float v)
  {
    cv_k = v / 20.; // = 0 to 5 range
  }

  static constexpr auto get_parameter (cv_tag)
  {
    return float_param ("%", 0.f, 100.f, 40.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<drive_tag, low_contour_tag, process_tag, cv_tag>;

private:
#if 0
  double get_slider_slider1()
  {
    // TODO: stub, add code for getting "slider2"
    // Range: min:0.0, max:10.0, default: 0.0, step: 0.1
    return 0.;
  }

  double get_slider_slider2()
  {
    // TODO: stub, add code for getting "slider2"
    // Range: min:0.0, max:10.0, default: 0.0, step: 0.1
    return 0.;
  }

  double get_slider_slider3()
  {
    // TODO: stub, add code for getting "slider3"
    // Range: min:0.0, max:10.0, default: 0.0, step: 0.1
    return 0.;
  }

  double get_slider_slider4()
  {
    // TODO: stub, add code for getting "slider4"
    // Range: min:0.0, max:10.0, default: 2.0, step: 0.1
    return 0.;
  }

  double get_slider_slider5()
  {
    // TODO: stub, add code for getting "slider5"
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.1
    return 0.;
  }

  double get_slider_slider6()
  {
    // TODO: stub, add code for getting "slider6"
    // Range: min:0.0, max:1.0, default: 1.0, step: None
    return 0.;
  }
#endif
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
  double a0_ap {};
  double a0_hp {};
  double a0_lp {};
  double a1_ap {};
  double a1_hp {};
  double a1_lp {};
  double a2_ap {};
  double a2_hp {};
  double a2_lp {};
  double alpha {};
  double b0_ap {};
  double b0_hp {};
  double b0_lp {};
  double b1_ap {};
  double b1_hp {};
  double b1_lp {};
  double b2_ap {};
  double b2_hp {};
  double b2_lp {};
  double cs {};
  double fc {};
  double fs {};
  double lp_fc {};
  double lp_r1pfc {};
  double nf_k {};
  double omega {};
  double pi2 {};
  double q {};
  double r_6 {};
  double rc0_k0 {};
  double rc0_k1 {};
  double sn {};
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double cv_k {};
  double g_hp {};
  double g_in {};
  double g_lp {};
  double g_out {};
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double buf_lp0 {};
  double buf_lp1 {};
  double cv_rc0 {};
  double cv_rc1 {};
  double out_rc0 {};
  double out_rc1 {};
  double x1_ap0 {};
  double x1_ap1 {};
  double x1_hp0 {};
  double x1_hp1 {};
  double x1_lp0 {};
  double x1_lp1 {};
  double x2_ap0 {};
  double x2_ap1 {};
  double x2_hp0 {};
  double x2_hp1 {};
  double x2_lp0 {};
  double x2_lp1 {};
  double y1_ap0 {};
  double y1_ap1 {};
  double y1_hp0 {};
  double y1_hp1 {};
  double y1_lp0 {};
  double y1_lp1 {};
  double y2_ap0 {};
  double y2_ap1 {};
  double y2_hp0 {};
  double y2_hp1 {};
  double y2_lp0 {};
  double y2_lp1 {};
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    g_in = g_out = 1.;
    set (drive_tag {}, 0.);
    set (low_contour_tag {}, 0.);
    set (process_tag {}, 0.);
    set (cv_tag {}, 40.);
    r_6      = 1. / 6.;
    nf_k     = std::pow (2., -96. * r_6);
    rc0_k0   = 0.0003;
    rc0_k1   = 1. - rc0_k0;
    lp_fc    = std::tan (3.141592653589793 * 0.48);
    lp_r1pfc = 1. / (1. + lp_fc);
    fc       = 700.;
    fs       = pc.get_sample_rate();
    pi2      = 2. * 3.141592653589793;
    q        = 0.23;
    omega    = pi2 * fc / fs;
    sn       = std::sin (omega);
    cs       = std::cos (omega);
    alpha    = sn / (2.0 * q);
    b0_ap    = 1.0 / (1.0 + alpha);
    b1_ap    = (-2.0 * cs) * b0_ap;
    b2_ap    = (1.0 - alpha) * b0_ap;
    a0_ap    = b2_ap;
    a1_ap    = b1_ap;
    a2_ap    = 1.;
    b0_lp    = 1.0 / (1.0 + alpha);
    b1_lp    = (-2.0 * cs) * b0_lp;
    b2_lp    = (1.0 - alpha) * b0_lp;
    a1_lp    = (1.0 - cs) * b0_lp;
    a0_lp    = a1_lp * 0.5;
    a2_lp    = a0_lp;
    b0_hp    = 1.0 / (1.0 + alpha);
    b1_hp    = (-2.0 * cs) * b0_hp;
    b2_hp    = (1.0 - alpha) * b0_hp;
    a1_hp    = -(1.0 + cs) * b0_hp;
    a0_hp    = -a1_hp * 0.5;
    a2_hp    = a0_hp;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    double cv0     = 0.;
    double cv1     = 0.;
    double in0     = 0.;
    double in1     = 0.;
    double nf0     = 0.;
    double nf1     = 0.;
    double out_ap0 = 0.;
    double out_ap1 = 0.;
    double out_dc0 = 0.;
    double out_dc1 = 0.;
    double out_hp0 = 0.;
    double out_hp1 = 0.;
    double out_lp0 = 0.;
    double out_lp1 = 0.;
    double out_sv0 = 0.;
    double out_sv1 = 0.;

    for (uint i = 0; i < samples; ++i) {
      auto& spl0 = outs[0][i];
      auto& spl1 = outs[1][i];
      spl0       = ins[0][i];
      spl1       = ins[1][i];
      in0        = spl0 * g_in;
      in1        = spl1 * g_in;
#if 1
      nf0 = (jsfx_rand (2.0) - 1.0) * nf_k;
      nf1 = (jsfx_rand (2.0) - 1.0) * nf_k;
#else
      nf0 = 0.;
      nf1 = 0.;
#endif
      out_ap0 = a0_ap * in0 + a1_ap * x1_ap0 + a2_ap * x2_ap0 - b1_ap * y1_ap0
        - b2_ap * y2_ap0;
      x2_ap0 = x1_ap0;
      x1_ap0 = in0;
      y2_ap0 = y1_ap0;
      y1_ap0 = out_ap0;

      out_ap1 = a0_ap * in1 + a1_ap * x1_ap1 + a2_ap * x2_ap1 - b1_ap * y1_ap1
        - b2_ap * y2_ap1;
      x2_ap1 = x1_ap1;
      x1_ap1 = in1;
      y2_ap1 = y1_ap1;
      y1_ap1 = out_ap1;

      out_lp0 = a0_lp * in0 + a1_lp * x1_lp0 + a2_lp * x2_lp0 - b1_lp * y1_lp0
        - b2_lp * y2_lp0;
      x2_lp0 = x1_lp0;
      x1_lp0 = in0;
      y2_lp0 = y1_lp0;
      y1_lp0 = out_lp0;

      out_lp1 = a0_lp * in1 + a1_lp * x1_lp1 + a2_lp * x2_lp1 - b1_lp * y1_lp1
        - b2_lp * y2_lp1;
      x2_lp1 = x1_lp1;
      x1_lp1 = in1;
      y2_lp1 = y1_lp1;
      y1_lp1 = out_lp1;

      out_hp0 = a0_hp * in0 + a1_hp * x1_hp0 + a2_hp * x2_hp0 - b1_hp * y1_hp0
        - b2_hp * y2_hp0;
      x2_hp0 = x1_hp0;
      x1_hp0 = in0;
      y2_hp0 = y1_hp0;
      y1_hp0 = out_hp0;

      out_hp1 = a0_hp * in1 + a1_hp * x1_hp1 + a2_hp * x2_hp1 - b1_hp * y1_hp1
        - b2_hp * y2_hp1;
      x2_hp1 = x1_hp1;
      x1_hp1 = in1;
      y2_hp1 = y1_hp1;
      y1_hp1 = out_hp1;

      cv0 = 1. - std::abs (in0) * cv_k;
      cv1 = 1. - std::abs (in1) * cv_k;

      cv_rc0 = rc0_k0 * cv0 + rc0_k1 * cv_rc0;
      cv_rc1 = rc0_k0 * cv1 + rc0_k1 * cv_rc1;

      out_sv0 = out_ap0 + out_hp0 * cv_rc0 * g_hp + out_lp0 * g_lp;
      out_sv1 = out_ap1 + out_hp1 * cv_rc1 * g_hp + out_lp1 * g_lp;

      out_rc0 = rc0_k0 * out_sv0 + rc0_k1 * out_rc0;
      out_rc1 = rc0_k0 * out_sv1 + rc0_k1 * out_rc1;
      out_dc0 = out_sv0 - out_rc0;
      out_dc1 = out_sv1 - out_rc1;

      out_lp0 = (buf_lp0 + lp_fc * out_dc0) * lp_r1pfc;
      buf_lp0 = lp_fc * (out_dc0 - out_lp0) + out_lp0;
      out_lp1 = (buf_lp1 + lp_fc * out_dc1) * lp_r1pfc;
      buf_lp1 = lp_fc * (out_dc1 - out_lp1) + out_lp1;

      spl0 = out_lp0 * g_out + nf0;
      spl1 = out_lp1 * g_out + nf1;
    }
  }
}; /* jsfx_process */

}} // namespace artv::liteon