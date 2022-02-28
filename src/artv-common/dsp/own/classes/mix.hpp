#pragma once

#include <array>
#include <cmath>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

namespace detail {
//----------------------------------------------------------------------------
static float to_normalized_balance (float v) // -100 to +100 to 0 to 1
{
  v = v / 100.;
  v = 0.5 * (v + 1.0);
  return v;
}
//------------------------------------------------------------------------------
static constexpr double smooth_time = 0.00025;
//------------------------------------------------------------------------------
struct gain_smoother {
public:
  //----------------------------------------------------------------------------
  void set_pan (float value)
  {
    value = detail::to_normalized_balance (value);
    if (pan == value) {
      return;
    }
    pan = value;
    update_target_gains();
  }
  //----------------------------------------------------------------------------
  void set_gain (float value)
  {
    if (gain == value) {
      return;
    }
    gain = value;
    update_target_gains();
  }
  //----------------------------------------------------------------------------
  void set_phaseinv_l (int value)
  {
    if (phase_inv_l == !!value) {
      return;
    }
    phase_inv_l = !!value;
    update_target_gains();
  }
  //----------------------------------------------------------------------------
  void set_phaseinv_r (int value)
  {
    if (phase_inv_r == !!value) {
      return;
    }
    phase_inv_r = !!value;
    update_target_gains();
  }
  //----------------------------------------------------------------------------
  void skip_smoothing()
  {
    l.setCurrentAndTargetValue (l.getTargetValue());
    r.setCurrentAndTargetValue (r.getTargetValue());
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    l.reset (pc.get_sample_rate(), smooth_time);
    r.reset (pc.get_sample_rate(), smooth_time);
  }
  //----------------------------------------------------------------------------
  bool is_smoothing() const { return l.isSmoothing() && r.isSmoothing(); }
  //----------------------------------------------------------------------------
  bool is_passthrough() const
  {
    return !is_smoothing() && gain == 1.f && pan == 0.5f && !phase_inv_l
      && !phase_inv_r;
  }
  //----------------------------------------------------------------------------
  juce::SmoothedValue<float> l {1.f};
  juce::SmoothedValue<float> r {1.f};
  //----------------------------------------------------------------------------
  float gain {1.f};
  float pan {.5f};
  bool  phase_inv_l {false};
  bool  phase_inv_r {false};
  //----------------------------------------------------------------------------
private:
  void update_target_gains() // -100 to 100
  {
    // see juce::dsp::Panner<T> p; // sin3dB
    float lval = std::sin (M_PI_2 * (1.0f - pan)) * M_SQRT2;
    float rval = std::sin (M_PI_2 * pan) * M_SQRT2;

    lval *= gain * (phase_inv_l ? -1.f : 1.f);
    rval *= gain * (phase_inv_r ? -1.f : 1.f);

    l.setTargetValue (lval);
    r.setTargetValue (rval);
  }
  //----------------------------------------------------------------------------
};
static constexpr auto gain_param = float_param ("", 0.f, 1.f, 1.0f, 0.001f);
static constexpr auto pan_param  = float_param ("", -100.f, 100.f, 0.0f, 0.01f);
static constexpr auto phinv_param
  = choice_param (0, make_cstr_array ("Off", "On"));
//------------------------------------------------------------------------------
struct ms_smoother {
public:
  void set_ms_ratio (float value)
  {
    value = detail::to_normalized_balance (value);
    if (ms == value) {
      return;
    }
    ms = value;
    update_target_ms();
  }
  //----------------------------------------------------------------------------
  void skip_smoothing()
  {
    m.setCurrentAndTargetValue (m.getTargetValue());
    s.setCurrentAndTargetValue (s.getTargetValue());
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    m.reset (pc.get_sample_rate(), smooth_time);
    s.reset (pc.get_sample_rate(), smooth_time);
  }
  //----------------------------------------------------------------------------
  bool is_smoothing() const { return m.isSmoothing() && s.isSmoothing(); }
  //----------------------------------------------------------------------------
  bool is_passthrough() const { return !is_smoothing() && ms == 0.5f; }
  //----------------------------------------------------------------------------
  juce::SmoothedValue<float> m {0.5f};
  juce::SmoothedValue<float> s {0.5f};
  float                      ms {.5f};

private:
  //----------------------------------------------------------------------------
  void update_target_ms() // -100 to 100
  {
    // linear decay on the stereo signal, as it usually is weak.
    auto  sval = std::min (ms, 0.5f);
    float mval;
    if (ms <= 0.5) {
      mval = 0.5f;
    }
    else {
      // exponential decay on the mono signal, as it usually is strong.
      // pow ((m - 1) * sqrt(2), 2)
      mval = ((ms - 1.) * M_SQRT2);
      mval *= mval;
    }
    m.setTargetValue (mval);
    s.setTargetValue (sval);
  }
};
//------------------------------------------------------------------------------
static constexpr auto ms_param = float_param ("", -100.f, 100.f, 0.f, 0.01f);
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
class dry_wet_mixer {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::mixer;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 2;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    for (auto& gs : _gain) {
      gs.reset (pc);
    }
    for (auto& mss : _ms) {
      mss.reset (pc);
    }
  }
  //----------------------------------------------------------------------------
  void skip_smoothing()
  {
    for (auto& gs : _gain) {
      gs.skip_smoothing();
    }
    for (auto& mss : _ms) {
      mss.skip_smoothing();
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block (crange<T*> outs, crange<T const*> ins, int samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= ((uint) bus_type));

    if (ins.size() == (uint) bus_type) {
      process_block_dry_only (outs, ins, samples);
    }
    else {
      process_block_dry_wet (outs, ins, samples);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_dry_wet (
    crange<T*>       outs,
    crange<T const*> ins,
    uint             samples)
  {
    using simd                  = juce::dsp::SIMDRegister<T>;
    constexpr size_t simd_elems = simd::SIMDNumElements;
    // TODO: test if/how the compiler vectorizes this and simplify if possible.
    // TODO: This is legacy. Use own vector wrappers instead.
    auto vect = [=] (std::array<T*, 4> c, uint blocks) {
      T* end = c[0] + (blocks * simd_elems);

      while (c[0] < end) {
        // dry channels
        auto dl  = simd::fromRawArray (c[0]);
        auto dr  = simd::fromRawArray (c[1]);
        auto reg = dl;
        dl       = dl + dr; // now L is M
        dr       = dr - reg; // now R is S
        dl *= _ms[dry].m.skip (simd_elems);
        dr *= _ms[dry].s.skip (simd_elems);
        reg = dl;
        dl  = dl - dr;
        dr  = reg + dr;
        dl *= _gain[dry].l.skip (simd_elems);
        dr *= _gain[dry].r.skip (simd_elems);
        // wet channels
        auto wl = simd::fromRawArray (c[2]);
        auto wr = simd::fromRawArray (c[3]);
        reg     = wl;
        wl      = wl + wr; // now L is M
        wr      = wr - reg; // now R is S
        wl *= _ms[wet].m.skip (simd_elems);
        wr *= _ms[wet].s.skip (simd_elems);
        reg = wl;
        wl  = wl - wr;
        wr  = reg + wr;
        wl *= _gain[wet].l.skip (simd_elems);
        wr *= _gain[wet].r.skip (simd_elems);
        // sum and global _gain
        dl += wl;
        dr += wr;
        dl *= _gain[global].l.skip (simd_elems);
        dr *= _gain[global].r.skip (simd_elems);

        dl.copyToRawArray (c[0]);
        dr.copyToRawArray (c[1]);
        c[0] += simd_elems;
        c[1] += simd_elems;
        c[2] += simd_elems;
        c[3] += simd_elems;
      }
    };

    auto unvect = [=] (std::array<T*, 4> c, uint rem) {
      T* end = c[0] + rem;
      while (c[0] < end) {
        // dry
        T dm = *c[0] + *c[1];
        T ds = *c[0] - *c[1];
        dm *= _ms[dry].m.getNextValue();
        ds *= _ms[dry].s.getNextValue();
        auto dl = (dm - ds) * _gain[dry].l.getNextValue();
        auto dr = (dm + ds) * _gain[dry].r.getNextValue();
        // wet
        T wm = *c[2] + *c[3];
        T ws = *c[2] - *c[3];
        wm *= _ms[wet].m.getNextValue();
        ws *= _ms[wet].s.getNextValue();
        auto wl = (wm - ws) * _gain[wet].l.getNextValue();
        auto wr = (wm + ws) * _gain[wet].r.getNextValue();
        // sum and global _gain
        *c[0] = (dl + wl) * _gain[global].l.getNextValue();
        *c[1] = (dr + wr) * _gain[global].r.getNextValue();

        ++c[0];
        ++c[1];
        ++c[2];
        ++c[3];
      }
    };
    for (uint i = 0; i < 2; ++i) {
      if (unlikely (ins[i] != outs[i])) {
        memcpy (outs[i], ins[i], samples * sizeof outs[i][0]);
      }
    }
    // this "const_cast" is because making "block_divide" const aware could be
    // a real mess. The wet channels are unnmodified.
    std::array<T*, 4> chnls
      = {outs[0], outs[1], const_cast<T*> (ins[2]), const_cast<T*> (ins[3])};

    block_divide (simd::SIMDRegisterSize, chnls, samples, vect, unvect);
  }
  //----------------------------------------------------------------------------
  // this can be invoked when you know there is no wet signal, it skips dry
  // panning, but it doesn't skip dry M/S balance.
  template <class T>
  void process_block_dry_only (
    crange<T*>       outs,
    crange<T const*> ins,
    uint             samples)
  {
    using simd                  = juce::dsp::SIMDRegister<T>;
    constexpr size_t simd_elems = simd::SIMDNumElements;

    if (
      _gain[global].is_passthrough() && _gain[dry].is_passthrough()
      && _ms[dry].is_passthrough()) {
      return;
    }

    // TODO: test if/how the compiler vectorizes this and simplify if
    // possible.
    // TODO: This is legacy. Use own vector wrappers instead.
    auto vect = [=] (std::array<T*, 2> c, uint blocks) {
      T* end = c[0] + (blocks * simd_elems);

      while (c[0] < end) {
        // dry channels
        auto dl  = simd::fromRawArray (c[0]);
        auto dr  = simd::fromRawArray (c[1]);
        auto reg = dl;
        dl       = dl + dr; // now L is M
        dr       = dr - reg; // now R is S
        dl *= _ms[dry].m.skip (simd_elems);
        dr *= _ms[dry].s.skip (simd_elems);
        reg = dl;
        dl  = dl - dr;
        dr  = reg + dr;
        dl *= _gain[global].l.skip (simd_elems);
        dr *= _gain[global].r.skip (simd_elems);

        dl.copyToRawArray (c[0]);
        dr.copyToRawArray (c[1]);
        c[0] += simd_elems;
        c[1] += simd_elems;
      }
    };

    auto unvect = [=] (std::array<T*, 2> c, uint rem) {
      T* end = c[0] + rem;
      while (c[0] < end) {
        // dry
        T dm = *c[0] + *c[1];
        T ds = *c[0] - *c[1];
        dm *= _ms[wet].m.getNextValue();
        ds *= _ms[wet].s.getNextValue();
        auto dl = (dm - ds) * _gain[global].l.getNextValue();
        auto dr = (dm + ds) * _gain[global].r.getNextValue();

        *c[0] = dl;
        *c[1] = dr;

        ++c[0];
        ++c[1];
      }
    };

    for (uint i = 0; i < 2; ++i) {
      if (unlikely (ins[i] != outs[i])) {
        memcpy (outs[i], ins[i], samples * sizeof outs[i][0]);
      }
    }
    // this "const_cast" is because making "block_divide" const aware could be
    // a real mess. The wet channels are unnmodified.
    std::array<T*, 2> dryb = {outs[0], outs[1]};
    block_divide (simd::SIMDRegisterSize, dryb, samples, vect, unvect);

    _gain[wet].l.skip (samples);
    _gain[wet].r.skip (samples);
    _gain[dry].l.skip (samples);
    _gain[dry].r.skip (samples);
    _ms[wet].m.skip (samples);
    _ms[wet].s.skip (samples);
  }
  //----------------------------------------------------------------------------
  struct pan_tag {};
  void set (pan_tag, float value) { _gain[global].set_pan (value); }
  static constexpr auto get_parameter (pan_tag) { return detail::pan_param; }
  //----------------------------------------------------------------------------
  struct gain_tag {};
  void set (gain_tag, float value) { _gain[global].set_gain (value); }
  static constexpr auto get_parameter (gain_tag) { return detail::gain_param; }
  //----------------------------------------------------------------------------
  struct phase_inv_l_tag {};
  void set (phase_inv_l_tag, int value)
  {
    _gain[global].set_phaseinv_l (value);
  }

  static constexpr auto get_parameter (phase_inv_l_tag)
  {
    return detail::phinv_param;
  }
  //----------------------------------------------------------------------------
  struct phase_inv_r_tag {};
  void set (phase_inv_r_tag, int value)
  {
    _gain[global].set_phaseinv_r (value);
  }

  static constexpr auto get_parameter (phase_inv_r_tag)
  {
    return detail::phinv_param;
  }
  //----------------------------------------------------------------------------
  struct dry_ms_ratio_tag {};
  void set (dry_ms_ratio_tag, float value) { _ms[dry].set_ms_ratio (value); }
  static constexpr auto get_parameter (dry_ms_ratio_tag)
  {
    return detail::ms_param;
  }
  //----------------------------------------------------------------------------
  struct dry_pan_tag {};
  void set (dry_pan_tag, float value) { _gain[dry].set_pan (value); }
  static constexpr auto get_parameter (dry_pan_tag)
  {
    return detail::pan_param;
  }
  //----------------------------------------------------------------------------
  struct dry_phase_inv_l_tag {};
  void set (dry_phase_inv_l_tag, int value)
  {
    _gain[dry].set_phaseinv_l (value);
  }

  static constexpr auto get_parameter (dry_phase_inv_l_tag)
  {
    return detail::phinv_param;
  }
  //----------------------------------------------------------------------------
  struct dry_phase_inv_r_tag {};
  void set (dry_phase_inv_r_tag, int value)
  {
    _gain[dry].set_phaseinv_r (value);
  }

  static constexpr auto get_parameter (dry_phase_inv_r_tag)
  {
    return detail::phinv_param;
  }
  //----------------------------------------------------------------------------
  struct wet_ms_ratio_tag {};
  void set (wet_ms_ratio_tag, float value) { _ms[wet].set_ms_ratio (value); }
  static constexpr auto get_parameter (wet_ms_ratio_tag)
  {
    return detail::ms_param;
  }
  //----------------------------------------------------------------------------
  struct wet_pan_tag {};
  void set (wet_pan_tag, float value) { _gain[wet].set_pan (value); }
  static constexpr auto get_parameter (wet_pan_tag)
  {
    return detail::pan_param;
  }
  //----------------------------------------------------------------------------
  struct wet_phase_inv_l_tag {};
  void set (wet_phase_inv_l_tag, int value)
  {
    _gain[wet].set_phaseinv_l (value);
  }

  static constexpr auto get_parameter (wet_phase_inv_l_tag)
  {
    return detail::phinv_param;
  }
  //----------------------------------------------------------------------------
  struct wet_phase_inv_r_tag {};
  void set (wet_phase_inv_r_tag, int value)
  {
    _gain[wet].set_phaseinv_r (value);
  }

  static constexpr auto get_parameter (wet_phase_inv_r_tag)
  {
    return detail::phinv_param;
  }
  //----------------------------------------------------------------------------
  struct dry_wet_ratio_tag {};
  void set (dry_wet_ratio_tag, float value)
  {
    value /= 100.;
    _gain[wet].set_gain (sqrt (value));
    _gain[dry].set_gain (sqrt (1. - value));
  }
  static constexpr auto get_parameter (dry_wet_ratio_tag)
  {
    return float_param ("", 0.f, 100.f, 100.f, 0.001f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    pan_tag,
    gain_tag,
    phase_inv_l_tag,
    phase_inv_r_tag,
    dry_ms_ratio_tag,
    dry_pan_tag,
    dry_phase_inv_l_tag,
    dry_phase_inv_r_tag,
    wet_ms_ratio_tag,
    wet_pan_tag,
    wet_phase_inv_l_tag,
    wet_phase_inv_r_tag,
    dry_wet_ratio_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  enum { dry, wet, global };
  std::array<detail::gain_smoother, 3> _gain;
  std::array<detail::ms_smoother, 2>   _ms;
};
//------------------------------------------------------------------------------
} // namespace artv
