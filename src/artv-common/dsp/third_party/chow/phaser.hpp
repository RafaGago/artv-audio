#pragma once

// Chow Phaser's DSP code is very-tightly coupled to JUCE, GUI objects, etc, so
// I unfortunately had to manually port it.
//
// From commit 9bf9cda60de10a591288c29a2ab2227cc971fe68
//
// It seems that both the original and the port don't behave well at high
// samplerates. Wild guess: Feedback delay changing? Does the feedback need a
// fractional delay line?

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace chow {
namespace detail {

//------------------------------------------------------------------------------
class chow_phaser_phase_section {
public:
  //----------------------------------------------------------------------------
  using value_type = double;
  //----------------------------------------------------------------------------
  static constexpr uint n_stages_max = 52;
  static constexpr uint order        = 1;
  //----------------------------------------------------------------------------
  enum coeffs { a1, b0, b1, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { n_states = n_stages_max };
  //----------------------------------------------------------------------------
  static void calc_coefs (xspan<value_type> c, float R, float fs)
  {
    assert (c.size() >= n_coeffs);
    // component values
    // R = jmax (1.0f, R);
    constexpr double C  = (float) 25e-9;
    double const     RC = R * C;

    // analog coefs
    double const b0s = RC;
    double const b1s = -1.0f;
    double const a0s = b0s;
    double const a1s = 1.0f;

    // bilinear transform
    double const K  = 2.0f * fs;
    double const a0 = a0s * K + a1s;
    c[b0]           = (b0s * K + b1s) / a0;
    c[b1]           = (-b0s * K + b1s) / a0;
    // a0: unused
    // c[a0]         = 1.0f;
    c[a1] = (-a0s * K + a1s) / a0;
  }
  //----------------------------------------------------------------------------
  static void reset_states (xspan<value_type> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static float tick (
    xspan<value_type const> c, // coeffs
    xspan<value_type>       s, // state
    value_type              x,
    float                   n_stages)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    // process integer stages
    for (uint stage = 0; stage < ((uint) n_stages); ++stage) {
      x = tick_stage (c, s, x, stage);
    }

    // process fractional stage
    float frac = n_stages - ((float) ((uint) n_stages));
    x = frac * tick_stage (c, s, x, (uint) n_stages) + (1.0f - frac) * x;

    return x;
  }
  //----------------------------------------------------------------------------
private:
  static value_type tick_stage (
    xspan<value_type const> c, // coeffs
    xspan<value_type>       s, // state,
    value_type              x,
    uint                    stage)
  {
    double y = s[stage] + x * c[b0];
    s[stage] = x * c[b1] - y * c[a1];

    return y;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
class chow_phaser_fb_section {
public:
  //----------------------------------------------------------------------------
  using value_type = double;
  //----------------------------------------------------------------------------
  enum coeffs {
    a1,
    a2,
    b0,
    b1,
    b2,
    n_own_coeffs,
    n_coeffs = n_own_coeffs + mystran_dc_blocker::n_coeffs
  };
  enum coeffs_int { n_coeffs_int };
  enum state {
    z1,
    z2,
    n_own_states,
    n_states = n_own_states + (3 * mystran_dc_blocker::n_states)
  };
  //----------------------------------------------------------------------------
  static void calc_coefs (
    xspan<value_type> c,
    double            R,
    double            fbAmt,
    double            fs)
  {
    using x1_vec = vec<value_type, 1>;

    assert (c.size() >= n_coeffs);
    // jassert (R > 0.0f);

    // component values
    // R = jmax (1.0f, R);
    constexpr double C  = (float) 15e-9;
    double const     RC = R * C;

    // analog coefs
    double const b0s = RC * RC;
    double const b1s = -2.0f * RC;
    double const b2s = 1.0f;
    double const a0s = b0s * (1.0f + fbAmt);
    double const a1s = -b1s * (1.0f - fbAmt);
    double const a2s = 1.0f + fbAmt;

    // frequency warping
    using fastmath   = juce::dsp::FastMathApproximations;
    double const wc  = calc_pole_freq (a0s, a1s, a2s);
    double const K   = (wc == 0.0f) ? 2.0f * fs : wc / tan (wc / (2.0 * fs));
    double const KSq = K * K;

    // bilinear transform
    double const a0 = a0s * KSq + a1s * K + a2s;
    // c[a0]          = a0 / a0;
    c[a1] = 2.0f * (a2s - a0s * KSq) / a0;
    c[a2] = (a0s * KSq - a1s * K + a2s) / a0;
    c[b0] = (b0s * KSq + b1s * K + b2s) / a0;
    c[b1] = 2.0f * (b2s - b0s * KSq) / a0;
    c[b2] = (b0s * KSq - b1s * K + b2s) / a0;

    // DC blockers, could be only done once...
    c.cut_head (n_own_coeffs);
    mystran_dc_blocker::reset_coeffs<x1_vec> (
      c.cast (x1_vec {}), make_vec (2.), 1. / fs);
  }
  //----------------------------------------------------------------------------
  static void reset_states (xspan<value_type> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static float tick (
    xspan<value_type const> c, // coeffs
    xspan<value_type>       s, // state
    value_type              x,
    float                   d1,
    float                   d2,
    float                   d3) noexcept
  {
    using x1_vec = vec<value_type, 1>;

    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    auto dc_coeffs = c.advanced (n_own_coeffs).cast<x1_vec const>();
    auto dc_states = s.advanced (n_own_states).cast (x1_vec {});

    double y      = s[z1] + x * c[b0];
    auto   ydrive = drive (y, d3);

    ydrive
      = mystran_dc_blocker::tick (dc_coeffs, dc_states, make_vec (ydrive))[0];

    s[z1] = drive (s[z2] + x * c[b1] - ydrive * c[a1], d1);
    dc_states.cut_head (mystran_dc_blocker::n_states);
    s[z1]
      = mystran_dc_blocker::tick (dc_coeffs, dc_states, make_vec (s[z1]))[0];

    s[z2] = drive (x * c[b2] - ydrive * c[a2], d2);
    dc_states.cut_head (mystran_dc_blocker::n_states);
    s[z2]
      = mystran_dc_blocker::tick (dc_coeffs, dc_states, make_vec (s[z2]))[0];
    return y;
  }
  //----------------------------------------------------------------------------
private:
  static double drive (double x, double drive) noexcept
  {
    return std::tanh (x * drive) / drive;
  }
  //----------------------------------------------------------------------------
  static double calc_pole_freq (double a, double b, double c)
  {
    auto radicand = b * b - 4. * a * c;
    if (radicand >= 0.0)
      return 0.0;

    return std::sqrt (-radicand) / (2.0 * a);
  }
};
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
class phaser {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::modulation;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct feedback_tag {};
  void                  set (feedback_tag, float v) { _params.fb = v; }
  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("", 0.0, 0.97, 0., 0.01);
  }
  //----------------------------------------------------------------------------
  struct modulation_tag {};
  void                  set (modulation_tag, float v) { _params.mod = v; }
  static constexpr auto get_parameter (modulation_tag)
  {
    return float_param ("", -1.0, 1.0, 0., 0.01);
  }
  //----------------------------------------------------------------------------
  struct lfo_freq_tag {};
  void set (lfo_freq_tag, float v)
  {

    v                = midi_note_to_hz (v); // FIX clang-format too many spaces
    _params.lfo_freq = v;
    reset_lfo();
  }

  static constexpr auto get_parameter (lfo_freq_tag)
  {
    return frequency_parameter_from_zero (16., 0.2);
  }
  //----------------------------------------------------------------------------
  struct lfo_depth_tag {};
  void                  set (lfo_depth_tag, float v) { _params.lfo_depth = v; }
  static constexpr auto get_parameter (lfo_depth_tag)
  {
    return float_param ("", 0.0, .95, 0.1, 0.01);
  }
  //----------------------------------------------------------------------------
  struct freq_mult_tag {};
  void set (freq_mult_tag, int v)
  {
    _params.freq_mult = (bool) v;
    reset_lfo();
  }

  static constexpr auto get_parameter (freq_mult_tag)
  {
    return choice_param (0, make_cstr_array ("1x", "10x"));
  }
  //----------------------------------------------------------------------------
  struct skew_tag {};
  void                  set (skew_tag, float v) { _params.skew = v; }
  static constexpr auto get_parameter (skew_tag)
  {
    return float_param ("", -3.0, 3.0, 0., 0.01);
  }
  //----------------------------------------------------------------------------
  struct stages_tag {};
  void                  set (stages_tag, float v) { _params.stages = v; }
  static constexpr auto get_parameter (stages_tag)
  {
    return float_param ("", 1., 50.0, 8., 0.01);
  }
  //----------------------------------------------------------------------------
  struct d1_tag {};
  void                  set (d1_tag, float v) { _params.d1 = v; }
  static constexpr auto get_parameter (d1_tag)
  {
    return float_param ("", 0.1, 5.0, 1., 0.01);
  }
  //----------------------------------------------------------------------------
  struct d2_tag {};
  void                  set (d2_tag, float v) { _params.d2 = v; }
  static constexpr auto get_parameter (d2_tag)
  {
    return float_param ("", 0.1, 5.0, 1., 0.01);
  }
  //----------------------------------------------------------------------------
  struct d3_tag {};
  void                  set (d3_tag, float v) { _params.d3 = v; }
  static constexpr auto get_parameter (d3_tag)
  {
    return float_param ("", 0.1, 5.0, 1., 0.01);
  }
  //----------------------------------------------------------------------------
  struct src_channel_tag {};
  void                  set (src_channel_tag, int v) { _params.src = v; }
  static constexpr auto get_parameter (src_channel_tag)
  {
    return choice_param (0, make_cstr_array ("Mid", "Left", "Right", "Side"));
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    feedback_tag,
    modulation_tag,
    lfo_freq_tag,
    lfo_depth_tag,
    freq_mult_tag,
    skew_tag,
    stages_tag,
    d1_tag,
    d2_tag,
    d3_tag,
    src_channel_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    using x1_t = vec<decltype (_smooth_coeff), 1>;

    _plugcontext = &pc;
    _t_spl       = 1.f / pc.get_sample_rate();

    smoother::reset_coeffs (
      xspan {&_smooth_coeff, 1}.cast (x1_t {}),
      vec_set<x1_t> (1. / 0.1),
      _t_spl);

    memset (&_ph_coeffs, 0, sizeof _ph_coeffs);
    memset (&_fb_coeffs, 0, sizeof _fb_coeffs);

    fb_section::reset_states (_fb_states);
    phase_section::reset_states (_ph_states);
    memset (&_param_state, 0, sizeof _param_state);
    _lfo.reset();

    _params.freq_mult = (bool) get_parameter (freq_mult_tag {}).defaultv;
    _params.lfo_depth = get_parameter (lfo_depth_tag {}).defaultv;
    _params.lfo_freq  = get_parameter (lfo_freq_tag {}).defaultv;
    _params.fb        = get_parameter (feedback_tag {}).defaultv;
    _params.mod       = get_parameter (modulation_tag {}).defaultv;
    _params.skew      = get_parameter (skew_tag {}).defaultv;
    _params.stages    = get_parameter (stages_tag {}).defaultv;
    _params.d1        = get_parameter (d1_tag {}).defaultv;
    _params.d2        = get_parameter (d2_tag {}).defaultv;
    _params.d3        = get_parameter (d3_tag {}).defaultv;
    _params.src       = get_parameter (src_channel_tag {}).defaultv;

    reset_lfo();
    param_targets_to_array (_param_state);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    alignas (sse_bytes) param_state_array param_tgt;
    alignas (sse_bytes) param_state_array smooth;
    param_targets_to_array (param_tgt);

    float sr = (float) _plugcontext->get_sample_rate();

    constexpr uint sse_step = vec_traits<f32_x4>().size;

    for (uint i = 0; i < samples; ++i) {
      // block LP parameter smoothing. 4 = SSE float
      for (uint j = 0; j < param_tgt.arr.size(); j += sse_step) {
        auto   spn  = xspan {&_param_state.arr[j], sse_step}.cast (f32_x4 {});
        auto   vecv = vec_load<f32_x4> (&param_tgt.arr[j]);
        f32_x4 out  = smoother::tick (xspan {&_smooth_coeff, 1}, spn, vecv);

        vec_store (&smooth.arr[j], out);
      }

      float lfo = _lfo.tick_sine()[0];
      lfo *= smooth.p.lfo_depth;
      float lfo_shape = light_shape (lfo, smooth.p.skew);

      constexpr float max_depth = 20.0f;
      auto            light_val = (max_depth + 0.1f) - (lfo_shape * max_depth);
      auto            rval = 100000.0f * std::pow (light_val / 0.1f, -0.75f);

      float mono;
      switch (_params.src) {
      case 0:
        mono = (ins[0][i] + ins[1][i]) * 0.5;
        break;
      case 1:
        mono = ins[0][i];
        break;
      case 2:
        mono = ins[1][i];
        break;
      case 3:
        mono = (ins[0][i] - ins[1][i]) * 0.5;
        break;
      default:
        break;
      }
      // when porting this I intuitively structured everything as I do with
      // filters: static functions + external arrays of coefficients and
      // states; for easy block-lowpass smoothing of coefficients. I didn't read
      // that this one was actually recalculating coefficients at sample rate.
      fb_section::calc_coefs (_fb_coeffs, rval, -1.f * smooth.p.fb, sr);
      double nomod_out = fb_section::tick (
        _fb_coeffs, _fb_states, mono, smooth.p.d1, smooth.p.d2, smooth.p.d3);

      phase_section::calc_coefs (_ph_coeffs, rval, sr);
      double mod_out = phase_section::tick (
        _ph_coeffs, _ph_states, nomod_out, smooth.p.stages);
      mod_out *= smooth.p.mod;
      mod_out += (1.0 - smooth.p.mod) * nomod_out;

      outs[0][i] = smooth.p.mod > 0.f ? mod_out : nomod_out;
      outs[1][i] = smooth.p.mod > 0.f ? nomod_out : mod_out;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  struct all_params {
    float fb;
    float mod;
    float lfo_freq;
    float lfo_depth;
    float skew;
    float stages;
    float d1;
    float d2;
    float d3;
    bool  freq_mult;
    uint  src;
  };

  struct smoothed_params {
    float fb;
    float mod;
    float lfo_depth;
    float skew;
    float stages;
    float d1;
    float d2;
    float d3;
  };

  union param_state_array {
    smoothed_params                                                         p;
    simd_array<float, sizeof (smoothed_params) / sizeof (float), sse_bytes> arr;
  };

  using smoother      = onepole_smoother;
  using phase_section = detail::chow_phaser_phase_section;
  using fb_section    = detail::chow_phaser_fb_section;
  //----------------------------------------------------------------------------
  void param_targets_to_array (param_state_array& pa)
  {
    pa.p.lfo_depth = _params.lfo_depth * (_params.freq_mult ? 0.8 : 1.);
    pa.p.d1        = _params.d1;
    pa.p.d2        = _params.d2;
    pa.p.d3        = _params.d3;
    pa.p.stages    = _params.stages;
    pa.p.skew      = std::pow (2.0f, _params.skew);
    pa.p.mod       = std::abs (_params.mod);
    pa.p.fb        = _params.fb;
  }
  //----------------------------------------------------------------------------
  void reset_lfo()
  {
    // the LFO coefficient is placed with the parameter smoothing....
    float mult = _params.freq_mult ? 10.f : 1.f;
    _lfo.set_freq (vec_set<1> (_params.lfo_freq * mult), _t_spl);
  }
  //----------------------------------------------------------------------------
  inline float light_shape (float x, float skewpow) noexcept
  {
    return (std::pow ((x + 1.0f) / 2.0f, skewpow) * 2.0f) - 1.0f;
  }
  //----------------------------------------------------------------------------
  std::array<double, phase_section::n_states> _ph_states;
  std::array<double, fb_section::n_states>    _fb_states;
  alignas (sse_bytes) param_state_array _param_state;

  std::array<double, phase_section::n_coeffs> _ph_coeffs;
  std::array<double, fb_section::n_coeffs>    _fb_coeffs;
  lfo<1>                                      _lfo;

  static_assert (smoother::n_states == 1, "");
  float           _smooth_coeff;
  all_params      _params;
  plugin_context* _plugcontext = nullptr;
  float           _t_spl;
};
//------------------------------------------------------------------------------
}} // namespace artv::chow
