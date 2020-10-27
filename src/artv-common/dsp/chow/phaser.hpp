#pragma once

// Chow Phaser's DSP code is very-tightly coupled to JUCE, GUI objects, etc, so
// I unfortunately had to manually port it.

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/dsp/own/blocks/filters/onepole.hpp"
#include "artv-common/dsp/own/blocks/oscillators/cheap_sine.hpp"
#include "artv-common/dsp/own/misc.hpp"
#include "artv-common/dsp/own/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace chow {
namespace detail {

//------------------------------------------------------------------------------
class chow_phaser_phase_section {
public:
  //----------------------------------------------------------------------------
  using value_type = float;
  //----------------------------------------------------------------------------
  static constexpr uint n_stages_max = 52;
  static constexpr uint order        = 1;
  //----------------------------------------------------------------------------
  enum coeffs { a1, b0, b1, n_coeffs };
  enum state { n_states = n_stages_max };
  //----------------------------------------------------------------------------
  static void calc_coefs (crange<value_type> c, float R, float fs)
  {
    assert (c.size() >= n_coeffs);
    // component values
    // R = jmax (1.0f, R);
    constexpr float C  = (float) 25e-9;
    const float     RC = R * C;

    // analog coefs
    const float b0s = RC;
    const float b1s = -1.0f;
    const float a0s = b0s;
    const float a1s = 1.0f;

    // bilinear transform
    const auto K  = 2.0f * fs;
    const auto a0 = a0s * K + a1s;
    c[b0]         = (b0s * K + b1s) / a0;
    c[b1]         = (-b0s * K + b1s) / a0;
    // a0: unused
    // c[a0]         = 1.0f;
    c[a1] = (-a0s * K + a1s) / a0;
  }
  //----------------------------------------------------------------------------
  static void reset_states (crange<value_type> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static float tick (
    crange<const value_type> c, // coeffs
    crange<value_type>       s, // state
    value_type               x,
    float                    n_stages)
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
    crange<const value_type> c, // coeffs
    crange<value_type>       s, // state,
    value_type               x,
    uint                     stage)
  {
    float y  = s[stage] + x * c[b0];
    s[stage] = x * c[b1] - y * c[a1];

    return y;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
class chow_phaser_fb_section {
public:
  //----------------------------------------------------------------------------
  using value_type = float;
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, b0, b1, b2, n_coeffs };
  enum state { z1, z2, n_states };
  //----------------------------------------------------------------------------
  static void calc_coefs (crange<value_type> c, float R, float fbAmt, float fs)
  {
    assert (c.size() >= n_coeffs);
    // jassert (R > 0.0f);

    // component values
    // R = jmax (1.0f, R);
    constexpr float C  = (float) 15e-9;
    const float     RC = R * C;

    // analog coefs
    const float b0s = RC * RC;
    const float b1s = -2.0f * RC;
    const float b2s = 1.0f;
    const float a0s = b0s * (1.0f + fbAmt);
    const float a1s = -b1s * (1.0f - fbAmt);
    const float a2s = 1.0f + fbAmt;

    // frequency warping
    using fastmath = juce::dsp::FastMathApproximations;
    const float wc = calc_pole_freq (a0s, a1s, a2s);
    const auto  K
      = (wc == 0.0f) ? 2.0f * fs : wc / fastmath::tan<float> (wc / (2.0f * fs));
    const auto KSq = K * K;

    // bilinear transform
    const float a0 = a0s * KSq + a1s * K + a2s;
    // c[a0]          = a0 / a0;
    c[a1] = 2.0f * (a2s - a0s * KSq) / a0;
    c[a2] = (a0s * KSq - a1s * K + a2s) / a0;
    c[b0] = (b0s * KSq + b1s * K + b2s) / a0;
    c[b1] = 2.0f * (b2s - b0s * KSq) / a0;
    c[b2] = (b0s * KSq - b1s * K + b2s) / a0;
  }
  //----------------------------------------------------------------------------
  static void reset_states (crange<value_type> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static float tick (
    crange<const value_type> c, // coeffs
    crange<value_type>       s, // state
    value_type               x,
    float                    d1,
    float                    d2,
    float                    d3) noexcept
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    float y      = s[z1] + x * c[b0];
    auto  ydrive = drive (y, d3);
    s[z1]        = drive (s[z2] + x * c[b1] - ydrive * c[a1], d1);
    s[z2]        = drive (x * c[b2] - ydrive * c[a2], d2);
    return y;
  }
  //----------------------------------------------------------------------------
private:
  static float drive (float x, float drive) noexcept
  {
    return std::tanh (x * drive) / drive;
  }
  //----------------------------------------------------------------------------
  static float calc_pole_freq (float a, float b, float c)
  {
    auto radicand = b * b - 4.0f * a * c;
    if (radicand >= 0.0f)
      return 0.0f;

    return std::sqrt (-radicand) / (2.0f * a);
  }
};
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
class phaser {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::modulation;
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
    reset_lfo_coeff();
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
    reset_lfo_coeff();
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
    _plugcontext = &pc;

    smoother::lowpass<float> (
      make_crange (_smooth_coeff), 1. / 0.1, pc.get_sample_rate());

    memset (&_ph_coeffs, 0, sizeof _ph_coeffs);
    memset (&_fb_coeffs, 0, sizeof _fb_coeffs);
    memset (&_lfo_coeffs, 0, sizeof _ph_coeffs);

    fb_section::reset_states (_fb_states);
    phase_section::reset_states (_ph_states);
    sin_osc::reset_states (_lfo_states);
    memset (&_param_state, 0, sizeof _param_state);

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

    reset_lfo_coeff();
    param_targets_to_array (_param_state);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    alignas (sse_bytes) param_state_array param_tgt;
    alignas (sse_bytes) param_state_array smooth;
    param_targets_to_array (param_tgt);

    float sr = (float) _plugcontext->get_sample_rate();

    constexpr uint sse_step = simd_flt::size;

    for (uint i = 0; i < samples; ++i) {
      // block LP parameter smoothing. 4 = SSE float
      for (uint j = 0; j < param_tgt.arr.size(); j += sse_step) {
        simd_flt out = smoother::tick_aligned<16, float> (
          make_crange (_smooth_coeff),
          make_crange (&_param_state.arr[j], sse_step),
          make_crange (&param_tgt.arr[j], sse_step));
        out.store_aligned (&smooth.arr[j]);
      }

      float lfo = sin_osc::tick (make_crange (smooth.p.lfo_coeff), _lfo_states);
      lfo *= smooth.p.lfo_depth;
      float lfo_shape = light_shape (lfo, smooth.p.skew);

      constexpr float max_depth = 20.0f;
      auto            light_val = (max_depth + 0.1f) - (lfo_shape * max_depth);
      auto            rval = 100000.0f * std::pow (light_val / 0.1f, -0.75f);

      float mono;
      switch (_params.src) {
      case 0:
        mono = (chnls[0][i] + chnls[1][i]) * 0.5;
        break;
      case 1:
        mono = chnls[0][i];
        break;
      case 2:
        mono = chnls[1][i];
        break;
      case 3:
        mono = (chnls[0][i] - chnls[1][i]) * 0.5;
        break;
      default:
        break;
      }
      // when porting this I intuitively structured everything as I do with
      // filters: static functions + external arrays of coefficients and
      // states; for easy block-lowpass smoothing of coefficients. I didn't read
      // that this one was actually recalculating coefficients at sample rate.
      fb_section::calc_coefs (_fb_coeffs, rval, -1.f * smooth.p.fb, sr);
      float nomod_out = fb_section::tick (
        _fb_coeffs, _fb_states, mono, smooth.p.d1, smooth.p.d2, smooth.p.d3);

      phase_section::calc_coefs (_ph_coeffs, rval, sr);
      float mod_out = phase_section::tick (
        _ph_coeffs, _ph_states, nomod_out, smooth.p.stages);
      mod_out *= smooth.p.mod;
      mod_out += (1.0f - smooth.p.mod) * nomod_out;

      chnls[0][i] = smooth.p.mod > 0.f ? mod_out : nomod_out;
      chnls[1][i] = smooth.p.mod > 0.f ? nomod_out : mod_out;
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
    float lfo_coeff;
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
  using sin_osc       = cheap_sin_osc;
  //----------------------------------------------------------------------------
  void param_targets_to_array (param_state_array& pa)
  {
    static_assert (sin_osc::n_coeffs == 1, "");

    pa.p.lfo_coeff = _lfo_coeffs[0];
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
  void reset_lfo_coeff()
  {
    // the LFO coefficient is placed with the parameter smoothing....
    float mult = _params.freq_mult ? 10.f : 1.f;
    sin_osc::calc_coefs (
      _lfo_coeffs, _params.lfo_freq * mult, _plugcontext->get_sample_rate());
  }
  //----------------------------------------------------------------------------
  inline float light_shape (float x, float skewpow) noexcept
  {
    return (std::pow ((x + 1.0f) / 2.0f, skewpow) * 2.0f) - 1.0f;
  }
  //----------------------------------------------------------------------------
  std::array<float, phase_section::n_states> _ph_states;
  std::array<float, fb_section::n_states>    _fb_states;
  std::array<float, sin_osc::n_states>       _lfo_states;
  alignas (sse_bytes) param_state_array _param_state;

  std::array<float, phase_section::n_coeffs> _ph_coeffs;
  std::array<float, fb_section::n_coeffs>    _fb_coeffs;
  std::array<float, sin_osc::n_coeffs>       _lfo_coeffs;

  static_assert (smoother::n_states == 1, "");
  float           _smooth_coeff;
  all_params      _params;
  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
}} // namespace artv::chow
