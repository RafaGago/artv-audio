#pragma once

// A very naive nonsense phaser. Still useful for phase rotation duties without
// modulation and feedback.

#include <array>
#include <cstdint>
#include <vector>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/biquad.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
class phaser {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::modulation;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct stages_tag {};
  void set (stages_tag, int v) { _pars.n_stages = v + 1; }

  static constexpr auto get_parameter (stages_tag)
  {
    // TODO: before release, not so many choices make a difference, reduce to
    // e.g. 2, 4, 6, 8, 12, 16, 24, 32
    return choice_param (
      6,
      make_cstr_array (
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16"),
      32);
  }
  //----------------------------------------------------------------------------
  struct stages_mode_tag {};
  void set (stages_mode_tag, int v)
  {
    bool lin          = v < m_total_modes;
    _pars.lin_lfo_mod = lin;
    _pars.mode        = v - (!lin * m_total_modes);
  }

  static constexpr auto get_parameter (stages_mode_tag)
  {
    return choice_param (
      8,
      make_cstr_array (
        "Lin Mod Exp F",
        "Lin Mod Exp F Spread",
        "Lin Mod Exp F Alternate",
        "Lin Mod Exp F Stereo",
        "Lin Mod Lin F",
        "Lin Mod Lin F Spread",
        "Lin Mod Lin F Alternate",
        "Lin Mod Lin F Stereo",
        "Exp Mod Exp F",
        "Exp Mod Exp F Spread",
        "Exp Mod Exp F Alternate",
        "Exp Mod Exp F Stereo",
        "Exp Mod Lin F",
        "Exp Mod Lin F Spread",
        "Exp Mod Lin F Alternate",
        "Exp Mod Lin F Stereo"),
      16);
  }

  enum stage_mode {
    m_exp,
    m_exp_spread,
    m_exp_alternate,
    m_exp_stereo,
    m_lin,
    m_lin_spread,
    m_lin_alternate,
    m_lin_stereo,
    m_total_modes,
  };

  static constexpr bool mode_is_exponential (uint v) { return v < m_lin; }
  //----------------------------------------------------------------------------
  struct lfo_rate_tag {};
  void set (lfo_rate_tag, float v)
  {
    v                 = midi_note_to_hz (v);
    _pars.lfo_hz_user = v;
    refresh_lfo_hz();
  }

  static constexpr auto get_parameter (lfo_rate_tag)
  {
    return frequency_parameter_from_zero (15.0, 0.66);
  }
  //----------------------------------------------------------------------------
  struct lfo_rate_sync_tag {};
  void set (lfo_rate_sync_tag, float v)
  {
    _pars.lfo_eights = v;
    refresh_lfo_hz();
  }

  static constexpr auto get_parameter (lfo_rate_sync_tag)
  {
    // NOTE: this has no 1/8 displayed because it breaks the display pane on
    // mix-matrix because of some (Juce?) bug.
    return float_param ("eights", 0., 32., 0., 1.);
  }
  //----------------------------------------------------------------------------
  struct lfo_depth_tag {};
  void set (lfo_depth_tag, float v)
  {
    v *= 0.01;
    _pars_smooth.target().lfo_depth = v;
  }

  static constexpr auto get_parameter (lfo_depth_tag)
  {
    return float_param ("%", 0., 100., 25., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lfo_wave_tag {};
  void set (lfo_wave_tag, int v) { _pars.lfo_wave = v; }

  static constexpr auto get_parameter (lfo_wave_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Sine", "Triangle", "S&H", "Noise", "Trapezoid", "Square", "Saw"),
      10);
  }
  //----------------------------------------------------------------------------
  struct lfo_start_phase_tag {};
  void set (lfo_start_phase_tag, float v)
  {
    _pars_smooth.target().lfo_start_phase = v;
  }

  static constexpr auto get_parameter (lfo_start_phase_tag)
  {
    return float_param ("deg", 0., 360., 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lfo_stereo_tag {};
  void set (lfo_stereo_tag, float v) { _pars_smooth.target().lfo_stereo = v; }

  static constexpr auto get_parameter (lfo_stereo_tag)
  {
    return float_param ("deg", 0., 360., 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct low_freq_tag {};
  void set (low_freq_tag, float v)
  {
    v                             = midi_note_to_hz (v);
    _pars_smooth.target().freq_lo = v;
  }

  static constexpr auto get_parameter (low_freq_tag)
  {
    return frequency_parameter (min_ap_freq, 10000., 225.);
  }
  //----------------------------------------------------------------------------
  struct high_freq_tag {};
  void set (high_freq_tag, float v)
  {
    v                             = midi_note_to_hz (v);
    _pars_smooth.target().freq_hi = v;
  }

  static constexpr auto get_parameter (high_freq_tag)
  {
    return frequency_parameter (min_ap_freq, 10000., 3300.);
  }
  //----------------------------------------------------------------------------
  static constexpr float max_feedback = 0.995;
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    bool neg = v < 0.f;
    v *= 0.01 * max_feedback;
    v = sqrt (abs (v));
    v = neg ? -v : v;

    _pars_smooth.target().feedback = v;
    scale_feedbacks (
      _pars_smooth.target().feedback, _pars_smooth.target().delay_feedback);
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("%", -100., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_hp_tag {};
  void set (feedback_hp_tag, float v)
  {
    _pars_smooth.target().feedback_hp = v * 0.01f;
  }

  static constexpr auto get_parameter (feedback_hp_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_lp_tag {};
  void set (feedback_lp_tag, float v)
  {
    _pars_smooth.target().feedback_lp = v * 0.01f;
  }

  static constexpr auto get_parameter (feedback_lp_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_sat_tag {};
  void set (feedback_sat_tag, float v)
  {
    v *= 0.01f;
    v = v * v * v;
    v *= 100.f;
    _pars_smooth.target().feedback_sat = v;
  }

  static constexpr auto get_parameter (feedback_sat_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct q_tag {};
  void set (q_tag, float v) { _pars_smooth.target().q = v; }

  static constexpr auto get_parameter (q_tag)
  {
    return float_param ("", 0.02, 3.5, 0.18, 0.001, 0.4);
  }
  //----------------------------------------------------------------------------
  struct parallel_mix_tag {};
  void set (parallel_mix_tag, float v)
  {
    _pars_smooth.target().parallel_mix = v * 0.005f;
  }

  static constexpr auto get_parameter (parallel_mix_tag)
  {
    return float_param ("%", -100., 100, 100., 0.1);
  }
  //----------------------------------------------------------------------------
  struct topology_tag {};
  void set (topology_tag, int v) { _pars.topology = v; }

  static constexpr auto get_parameter (topology_tag)
  {
    return choice_param (1, make_cstr_array ("1 pole", "2 pole"), 16);
  }

  enum topologies { t_1_pole, t_2_pole };
  //----------------------------------------------------------------------------
  struct delay_feedback_tag {};
  void set (delay_feedback_tag, float v)
  {
    bool neg = v < 0.f;
    v *= 0.01 * 0.96;
    v = sqrt (abs (v));
    v = neg ? -v : v;

    _pars_smooth.target().delay_feedback = v;
    scale_feedbacks (
      _pars_smooth.target().feedback, _pars_smooth.target().delay_feedback);
  }

  static constexpr auto get_parameter (delay_feedback_tag)
  {
    return float_param ("%", -100., 100, 0., 0.01);
  }
  //----------------------------------------------------------------------------
  static constexpr float max_delay_sec = 0.0015f;

  struct delay_time_tag {};
  void set (delay_time_tag, float v)
  {
    v *= 0.01f;
    v *= max_delay_sec * _srate;
    _pars_smooth.target().delay_spls = v;
  }

  static constexpr auto get_parameter (delay_time_tag)
  {
    return float_param ("%", 0., 100, 33., 0.001);
  }
  //----------------------------------------------------------------------------
  static constexpr float max_delay_factor    = 4.f;
  static constexpr float ln_max_delay_factor = M_LN2 * 2.f;

  struct delay_lfo_tag {};
  void set (delay_lfo_tag, float v)
  {
    v *= 0.01f;
    _pars_smooth.target().delay_lfo = v;
  }

  static constexpr auto get_parameter (delay_lfo_tag)
  {
    return float_param ("%", -100., 100, 0., 0.001);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    _srate       = pc.get_sample_rate();
    _t_spl       = 1.f / _srate;
    _allpass1p.reset_states_cascade();
    _allpass2p.reset_states_cascade();
    _feedback_shelf.reset_states_cascade();
    _feedback_delay_shelf.reset_states_cascade();
    _dc_blocker.reset_states();
    _lfos.reset();
    _pars_smooth.reset (_t_spl, 1.f / 0.01f);

    // ensure that all parameter's new values will be different
    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    _pars_smooth.set_all_from_target();

    refresh_lfo_hz();

    _dc_blocker.reset_coeffs (vec_set<4> (2.f), _t_spl);

    _n_processed_samples = 0;

    // Sample rates 44100 multiples update every 362.811us
    uint sr_order      = get_samplerate_order (pc.get_sample_rate()) + 3;
    _control_rate_mask = lsb_mask<uint> (sr_order);
    _feedback_samples  = vec_set<f64_x2> (0.);
    _lfo_t_spl         = _t_spl * (_control_rate_mask + 1);

    // setting up delay
    reset_memory_parts();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {

      constexpr uint vec_size = vec_traits<f32_x4>().size;

      // parameter smoothing (could be done at a subblock/control rate).
      _pars_smooth.tick();

      all_parameters pars;
      static_assert (
        onepole_smoother::n_states == 1,
        "For the assignment below to work only the Z-1 state has to exist");
      *((smoothed_parameters*) &pars)   = _pars_smooth.get();
      *((unsmoothed_parameters*) &pars) = _pars;

      // control rate/mdoulation block -----------------------------------------
      if ((_n_processed_samples & _control_rate_mask) == 0) {
        run_control_rate (pars);
      }
      // sample-wise processing block (unbreakable, it has 1 sample delays) ----
      f32_x4 fb = {
        (float) _feedback_samples[0],
        (float) _feedback_samples[1],
        (float) _feedback_samples[0],
        (float) _feedback_samples[1]};
      // regular processing
      f32_x4 out {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      out += fb * pars.feedback;
      f64_x2 feedforward {out[0], out[1]};

      auto n_spls = pars.delay_spls;
      n_spls *= exp (ln_max_delay_factor * pars.lfo_last * pars.delay_lfo);
      auto delayed = _delay.get (n_spls, 0);
      out += delayed * pars.delay_feedback;

      if (pars.topology == t_2_pole) {
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass2p.tick_on_idx (g, out);
        }
      }
      else {
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass1p.tick_on_idx (g, out);
        }
      }
      out         = _dc_blocker.tick (out);
      auto del_fb = _feedback_delay_shelf.tick_on_idx (k_shelf_lo, out);
      del_fb      = _feedback_delay_shelf.tick_on_idx (k_shelf_hi, out);
      _delay.push (xspan {&del_fb, 1});

      f64_x2 outx2    = {(double) out[0], (double) out[1]};
      f64_x2 parallel = {(double) out[2], (double) out[3]};

      outx2 *= 0.5f + (0.5f - abs (pars.parallel_mix));
      outx2 += pars.parallel_mix * parallel;

      auto main_fb = outx2;
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_lo, main_fb);
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_hi, main_fb);
      main_fb = main_fb / vec_sqrt (main_fb * main_fb * pars.feedback_sat + 1.);
      _feedback_samples = main_fb;

      outs[0][i] = outx2[0];
      outs[1][i] = outx2[1];
    }
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    stages_tag,
    stages_mode_tag,
    lfo_rate_tag,
    lfo_depth_tag,
    lfo_start_phase_tag,
    lfo_stereo_tag,
    lfo_wave_tag,
    low_freq_tag,
    high_freq_tag,
    feedback_tag,
    q_tag,
    parallel_mix_tag,
    delay_feedback_tag,
    feedback_hp_tag,
    feedback_lp_tag,
    feedback_sat_tag,
    delay_time_tag,
    delay_lfo_tag>;
  //----------------------------------------------------------------------------
private:
  static constexpr uint  max_stages    = 16;
  static constexpr float min_ap_freq   = 20.f;
  static constexpr uint  n_channels    = 2;
  static constexpr uint  n_ap_channels = 4;
  //----------------------------------------------------------------------------
  // As feedbacks are smoothed, this might need to be done sample-wise.
  void scale_feedbacks (float& g1, float& g2)
  {
    // keep the sum of both feebacks below unity
    float fb_gain = abs (g1) + abs (g2);
    float att     = (fb_gain <= max_feedback) ? 1.f : (max_feedback / fb_gain);
    g1 *= att;
    g2 *= att;
  }
  //----------------------------------------------------------------------------
  void reset_memory_parts()
  {
    auto delay_size
      = pow2_round_ceil ((uint) (_srate * max_delay_sec * max_delay_factor));
    delay_size = _delay.n_required_elems (delay_size, 1);
    delay_size = div_ceil<unsigned> (
      delay_size, sizeof _delay_mem[0] / sizeof _delay.get (0.f, 0));

    _delay_mem.clear();
    _delay_mem.resize (delay_size);
    auto mem = xspan {_delay_mem};
    _delay.reset (mem, 1);
  }
  //----------------------------------------------------------------------------
  void refresh_lfo_hz()
  {
    float hz = 0.f;
    if (_pars.lfo_eights != 0.f) {
      hz = (8. * _plugcontext->get_play_state().bpm) / (60. * _pars.lfo_eights);
    }
    hz += _pars.lfo_hz_user;
    _pars_smooth.target().lfo_hz_final = hz;
  }
  //----------------------------------------------------------------------------
  template <class T>
  T freq_to_delay_spls (T freq)
  {
    if constexpr (is_vec_v<T>) {
      using VT = vec_value_type_t<T>;
      return (VT) _srate / ((VT) 2 * freq);
    }
    else {
      return (T) _srate / ((T) 2 * freq);
    }
  }
  //----------------------------------------------------------------------------
  struct all_parameters;
  void run_control_rate (all_parameters const& pars)
  {
    // this is very naive, note that we could e.g. be interpolating SVF
    // coefficients internally.
    constexpr uint vec_size = vec_traits<f32_x4>().size;

    _lfos.set_freq (vec_set<2> (pars.lfo_hz_final), _lfo_t_spl);
    auto stereo_ph
      = phase<1> {phase_tag::degrees {}, pars.lfo_stereo}.get_raw (0);
    if (pars.lfo_hz_final != 0.f) {
      // updating stereo phase diff.
      auto ph = _lfos.get_phase();
      ph.set_raw (ph.get_raw (0) + stereo_ph, 1);
      _lfos.set_phase (ph);
    }
    else {
      auto start_ph = phase<n_channels> {
        phase_tag::degrees {}, vec_set<n_channels> (pars.lfo_start_phase)};
      start_ph.set_raw (start_ph.get_raw (0) + stereo_ph, 1);
      _lfos.set_phase (start_ph);
    }
    vec<float, 2> lfov;
    switch (pars.lfo_wave) {
    case 0:
      lfov = _lfos.tick_sine();
      break;
    case 1:
      lfov = _lfos.tick_triangle();
      break;
    case 2:
      lfov = _lfos.tick_sample_hold();
      break;
    case 3:
      lfov = _lfos.tick_filt_sample_and_hold();
      break;
    case 4:
      lfov = _lfos.tick_trapezoid (vec_set<2> (0.3f));
      break;
    case 5:
      lfov = _lfos.tick_square();
      break;
    case 6:
      lfov = _lfos.tick_saw();
      break;
    }
    _pars_smooth.target().lfo_last      = lfov[0];
    auto                  interp_stages = (double) (pars.n_stages * 2) - 1;
    std::array<double, 2> fconstant;

    auto freq_hi = std::max (pars.freq_hi, pars.freq_lo);
    auto freq_lo = std::min (pars.freq_hi, pars.freq_lo);

    switch (pars.mode) {
    case m_exp:
    case m_exp_spread:
    case m_exp_alternate:
      fconstant[0] = 1.;
      fconstant[1] = pow (freq_hi / freq_lo, 1. / interp_stages);
      break;
    case m_exp_stereo:
      interp_stages *= 2;
      fconstant[0] = pow (freq_hi / freq_lo, 1. / interp_stages);
      fconstant[1] = fconstant[0];
      break;
    case m_lin:
    case m_lin_spread:
    case m_lin_alternate:
      fconstant[0] = 0.f;
      fconstant[1] = (freq_hi - freq_lo) / interp_stages;
      break;
    case m_lin_stereo:
      interp_stages *= 2;
      fconstant[0] = (freq_hi - freq_lo) / interp_stages;
      fconstant[1] = fconstant[0];
      break;
    default:
      fconstant[0] = 0.;
      fconstant[1] = 0.;
      break;
    }
    auto   f         = freq_lo;
    auto   lfo_depth = pars.lfo_depth;
    f32_x4 q_fact {lfov[0], lfov[1], lfov[1], lfov[0]};
    q_fact = vec_exp (q_fact * -lfo_depth * 0.23f);

    for (uint s = 0; s < pars.n_stages; ++s) {
      f32_x4 freqs {};
      f32_x4 qs = vec_set<f32_x4> (pars.q);
      // fill freqs
      for (uint ap = 0; ap < (vec_size / 2); ++ap) {
        float max_depth = pars.lin_lfo_mod ? 0.6 : 1.7;
        float k;
        uint  stage = (s * 4) + ap;

        switch (pars.mode) {
        case m_exp:
        case m_lin:
        case m_exp_stereo:
        case m_lin_stereo:
          k = max_depth;
          break;
        case m_exp_spread:
        case m_lin_spread:
          k = (stage < (pars.n_stages * 2)) ? max_depth : -max_depth;
          break;
        case m_exp_alternate:
        case m_lin_alternate:
          k = ((stage / 2) & 1) ? max_depth : -max_depth;
          break;
        default:
          break;
        }

        std::array<decltype (f), 2> fs;
        for (uint iv = 0; iv < fs.size(); ++iv) {
          fs[iv] = f;
          if (mode_is_exponential (pars.mode)) {
            f *= fconstant[iv];
          }
          else {
            f += fconstant[iv];
          }
        }

        if (pars.lin_lfo_mod) {
          freqs[ap * 2]     = fs[0] + (lfov[0] * lfo_depth * fs[0] * k);
          freqs[ap * 2 + 1] = fs[1] + (lfov[1] * lfo_depth * fs[1] * k);
        }
        else {
          freqs[ap * 2]     = fs[0] * exp (lfov[0] * lfo_depth * k);
          freqs[ap * 2 + 1] = fs[1] * exp (lfov[1] * lfo_depth * k);
        }
      }
      freqs
        = vec_clamp (freqs, min_ap_freq, std::min (_srate * 0.4999f, 22000.f));
      if (pars.topology == t_1_pole) {
        _allpass1p.reset_coeffs_on_idx (s, freqs, _t_spl);
      }
      else if (pars.topology == t_2_pole) {
        _allpass2p.reset_coeffs_on_idx (s, freqs, qs, _t_spl);
      }
      if (s == 0) {
        constexpr float cut_ratio = 0.3f;
        constexpr float gain_db   = -20.f;

        auto cutfreq = f64_x2 {freqs[0], freqs[1]};
        cutfreq -= pars.feedback_hp * cutfreq * cut_ratio;
        _feedback_shelf.reset_coeffs_on_idx (
          k_shelf_lo,
          vec_max (210., cutfreq),
          vec_set<f64_x2> (0.35),
          vec_set<f64_x2> (pars.feedback_hp * gain_db),
          _t_spl,
          lowshelf_tag {});

        auto del_cutfreq = freqs;
        del_cutfreq -= pars.feedback_hp * freqs * cut_ratio;
        _feedback_delay_shelf.reset_coeffs_on_idx (
          k_shelf_lo,
          vec_max (210.f, del_cutfreq),
          vec_set<f32_x4> (0.35),
          vec_set<f32_x4> (pars.feedback_hp * gain_db),
          _t_spl,
          lowshelf_tag {});
      }
      if (s == (pars.n_stages - 1)) {
        constexpr float cut_ratio = 0.3f;
        constexpr float gain_db   = -6.f;

        auto cutfreq = f64_x2 {freqs[0], freqs[1]};
        cutfreq += pars.feedback_lp * cutfreq * cut_ratio;
        _feedback_shelf.reset_coeffs_on_idx (
          k_shelf_hi,
          vec_max (4000., cutfreq),
          vec_set<f64_x2> (0.35),
          vec_set<f64_x2> (pars.feedback_lp * gain_db),
          _t_spl,
          highshelf_tag {});

        auto del_cutfreq = freqs;
        del_cutfreq -= pars.feedback_hp * freqs * cut_ratio;
        _feedback_delay_shelf.reset_coeffs_on_idx (
          k_shelf_hi,
          vec_max (4000.f, del_cutfreq),
          vec_set<f32_x4> (0.35),
          vec_set<f32_x4> (pars.feedback_lp * gain_db),
          _t_spl,
          highshelf_tag {});
      }
    }
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    float lfo_eights;
    float lfo_hz_user;
    uint  lfo_wave;
    uint  n_stages;
    uint  topology;
    uint  mode;
    bool  lin_lfo_mod;
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float freq_lo;
    float freq_hi;
    float q;
    float lfo_hz_final;

    float lfo_stereo;
    float lfo_depth;
    float lfo_start_phase;
    float feedback;

    float delay_spls;
    float delay_feedback;
    float delay_lfo;
    float lfo_last;

    float parallel_mix;
    float feedback_hp;
    float feedback_lp;
    float feedback_sat;
  };
  //----------------------------------------------------------------------------
  struct all_parameters : public unsmoothed_parameters,
                          public smoothed_parameters {};
  //----------------------------------------------------------------------------
  f64_x2                _feedback_samples;
  unsmoothed_parameters _pars;

  std::vector<f32_x4, overaligned_allocator<f32_x4, 128>> _delay_mem;

  value_smoother<float, smoothed_parameters> _pars_smooth;
  // TODO: use a variant?
  part_class_array<andy::svf_allpass, f32_x4, max_stages> _allpass2p;
  part_class_array<onepole_allpass, f32_x4, max_stages>   _allpass1p;
  enum { k_shelf_lo, k_shelf_hi };
  part_class_array<andy::svf, f64_x2, 2> _feedback_shelf;
  part_class_array<andy::svf, f32_x4, 2> _feedback_delay_shelf;

  modulable_thiran1_delay_line<f32_x4, 1>      _delay {};
  part_class_array<mystran_dc_blocker, f32_x4> _dc_blocker;

  lfo<n_channels> _lfos; // 0 = L, 1 = R
  uint            _n_processed_samples;
  uint            _control_rate_mask;
  float           _lfo_t_spl;
  float           _t_spl;
  float           _srate;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
