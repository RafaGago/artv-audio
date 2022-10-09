#pragma once

#include <algorithm>
#include <gcem.hpp>
#include <numeric>
#include <vector>

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/diffusion_matrix.hpp"
#include "artv-common/dsp/own/classes/ducker.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/reverb_tools.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/filters/saike.hpp"
#include "artv-common/dsp/own/parts/interpolation/stateless.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/third_party/saike/transience.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

// Both thiran and sinc of different sizes consume the same CPU, the sound
// is not that bad on the Thiran ones and they don't require tables, so favoring
// those.

#define DIFFUSE_DELAY_USE_THIRAN_TAPS      1
#define DIFFUSE_DELAY_USE_SINC_TAPS        0
#define DIFFUSE_DELAY_USE_THIRAN_DIFFUSORS 1
namespace artv {

//------------------------------------------------------------------------------
class diffuse_delay {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::delay;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  static constexpr uint max_t_beats = 2;
  //----------------------------------------------------------------------------
  struct sixteenths_tag {};
  void set (sixteenths_tag, float v)
  {
    v = v * (float) (1 / 16.) * _param.spls_x_beat;
    if (v == _extpar.delay_spls) {
      return;
    }
    _extpar.delay_spls = v;
    update_taps();
  }

  static constexpr auto get_parameter (sixteenths_tag)
  {
    return float_param ("sixteenths", 0.5f, max_t_beats * 16, 6.f, 0.01f, 0.5f);
  }
  //----------------------------------------------------------------------------
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    v = v * 0.01f;
    v *= v;
    v = 0.004f + v * 0.996f;
    if (v == _extpar.feedback) {
      return;
    }
    _extpar.feedback = v;
    update_taps();
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("%", 0.f, 100.f, 30.f, 0.001f);
  }
  //----------------------------------------------------------------------------
  struct diffusion_tag {};
  void set (diffusion_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.diffusion) {
      return;
    }
    _extpar.diffusion      = v;
    float weight           = cosf (v * (float) (M_PI / 4.));
    _param.mtx_angle[0][0] = weight;
    _param.mtx_angle[0][1] = sqrt (1.f - weight * weight);
  }

  static constexpr auto get_parameter (diffusion_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct ducking_threshold_tag {};
  void set (ducking_threshold_tag, float v)
  {
    if (v == _extpar.ducking_threshold) {
      return;
    }
    _extpar.ducking_threshold = v;
    _ducker.set_threshold (vec_set<2> ((double) v));
  }

  static constexpr auto get_parameter (ducking_threshold_tag)
  {
    return float_param ("dB", -40.f, 12.f, 12.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct ducking_speed_tag {};
  void set (ducking_speed_tag, float v)
  {
    if (v == _extpar.ducking_speed) {
      return;
    }
    _extpar.ducking_speed = v;
    v *= 0.01f;
    _ducker.set_speed (vec_set<f64_x2> (v * v), t_spl);
  }

  static constexpr auto get_parameter (ducking_speed_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct gain_tag {};
  void set (gain_tag, float v)
  {
    if (v == _extpar.gain) {
      return;
    }
    _extpar.gain     = v;
    _param.main_gain = db_to_gain (v);
  }

  static constexpr auto get_parameter (gain_tag)
  {
    return float_param ("dB", -50.f, 12.f, -10.f, 0.2f);
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void set (mode_tag, int v)
  {
    if (v == _extpar.mode) {
      return;
    }
    _extpar.mode = v;
  }

  enum {
    m_stereo,
    m_ping_pong,
    m_ping_pong_stereo,
    m_123,
    m_132,
    m_1223,
    m_1234,
    m_2413,
    m_1324,
    m_2314,
    m_1423,
    m_chorus
  };

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      1,
      make_cstr_array (
        "Stereo",
        "Ping-Pong",
        "Ping-Pong Stereo",
        "L-C-R",
        "L-R-C",
        "L-C-C-R",
        "L-LC-RC-R",
        "LC-R-L-RC",
        "L-RC-LC-R",
        "LC-RC-L-R",
        "L-R-LC-RC",
        "4 synced taps"),
      32);
  }
  //----------------------------------------------------------------------------
  struct mod_freq_tag {};
  void set (mod_freq_tag, float v)
  {
    if (v == _extpar.mod_freq) {
      return;
    }
    _extpar.mod_freq = v;
    _mod_lfo.set_freq (vec_set<n_taps> (v), t_spl);
  }

  static constexpr auto get_parameter (mod_freq_tag)
  {
    return float_param ("Hz", 0.f, 12.f, 0.17f, 0.001f, 0.4f);
  }
  //----------------------------------------------------------------------------
  struct mod_depth_tag {};
  void set (mod_depth_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.mod_depth) {
      return;
    }
    _extpar.mod_depth = v;
  }

  static constexpr auto get_parameter (mod_depth_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct mod_mode_tag {};
  void set (mod_mode_tag, uint v)
  {
    if (v == _extpar.mod_mode) {
      return;
    }
    _extpar.mod_mode = v;
  }

  static constexpr auto get_parameter (mod_mode_tag)
  {
    return choice_param (
      0, make_cstr_array ("Random", "Sine", "Triangle", "Trapezoid"), 10);
  }
  //----------------------------------------------------------------------------
  struct damp_tag {};
  void set (damp_tag, float v)
  {
    v *= 0.01f;
    v = 1.f - v;
    if (v == _extpar.damp_ratio) {
      return;
    }
    _extpar.damp_ratio = v;
    update_damp();
  }

  static constexpr auto get_parameter (damp_tag)
  {
    return float_param ("%", 0.f, 100.f, 15.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct freq_spread_tag {};
  void set (freq_spread_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.freq_spread) {
      return;
    }
    _extpar.freq_spread = v;
    update_damp();
    update_bp();
  }

  static constexpr auto get_parameter (freq_spread_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct tilt_db_tag {};
  void set (tilt_db_tag, float v)
  {
    if (v == _extpar.tilt_db) {
      return;
    }
    _extpar.tilt_db = v;
    _tilt.reset_coeffs (
      vec_set<2> (330.),
      vec_set<2> (0.5),
      vec_set<2> ((double) -v),
      (double) t_spl);
  }

  static constexpr auto get_parameter (tilt_db_tag)
  {
    return float_param ("dB", -14.f, 14.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct desync_tag {};
  void set (desync_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.desync) {
      return;
    }
    _extpar.desync = v;
    update_taps();
  }

  static constexpr auto get_parameter (desync_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct transients_tag {};
  void set (transients_tag, float v)
  {
    v *= 0.008f;
    if (v == _extpar.transients) {
      return;
    }
    _extpar.transients = v;
    _transients.set (saike::transience::strength_tag {}, v);
    _transients.set (saike::transience::strength2_tag {}, -v);
  }

  static constexpr auto get_parameter (transients_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct hipass_tag {};
  void set (hipass_tag, float v)
  {
    v *= 0.01f;
    v *= v;
    if (v == _extpar.hp) {
      return;
    }
    _extpar.hp = v;
    _filters.reset_coeffs<hp_idx> (vec_set<n_taps> (4.f + v * 396.f), t_spl);
  }

  static constexpr auto get_parameter (hipass_tag)
  {
    return float_param ("%", 0.f, 100.f, 25.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct bp_drive_tag {};
  void set (bp_drive_tag, float v)
  {
    constexpr float dbrange = 11.f;
    v *= (0.01 * dbrange);
    if (v == _extpar.bp_drive_db) {
      return;
    }
    _extpar.bp_drive_db = v;
    _extpar.bp_drive    = db_to_gain (v);
  }

  static constexpr auto get_parameter (bp_drive_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct bp_freq_tag {};
  void set (bp_freq_tag, float v)
  {
    if (v == _extpar.bp_note) {
      return;
    }
    _extpar.bp_note = v;
    update_bp();
  }

  static constexpr float bp_min_hz = 100.;
  static constexpr float bp_max_hz = 3500.;

  static constexpr auto get_parameter (bp_freq_tag)
  {
    return frequency_parameter (bp_min_hz, bp_max_hz, 440.0);
  }
  //----------------------------------------------------------------------------
  struct bp_envfollow_tag {};
  void set (bp_envfollow_tag, float v) { _extpar.bp_envfollow = v * 0.01; }

  static constexpr auto get_parameter (bp_envfollow_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct bp_wet_dry_tag {};
  void set (bp_wet_dry_tag, float v)
  {
    v *= 0.01;
    if (v == _extpar.bp_wetdry) {
      return;
    }
    _extpar.bp_wetdry = v;
    update_bp();
  }

  static constexpr auto get_parameter (bp_wet_dry_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct bp_reso_tag {};
  void set (bp_reso_tag, float v)
  {
    v *= 0.01;
    _extpar.bp_reso = 0.0001f + v * 0.9999f;
  }

  static constexpr auto get_parameter (bp_reso_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.001f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    feedback_tag,
    gain_tag,
    sixteenths_tag,
    damp_tag,
    freq_spread_tag,
    tilt_db_tag,
    mode_tag,
    diffusion_tag,
    desync_tag,
    mod_mode_tag,
    mod_freq_tag,
    mod_depth_tag,
    ducking_speed_tag,
    ducking_threshold_tag,
    transients_tag,
    hipass_tag,
    bp_drive_tag,
    bp_freq_tag,
    bp_reso_tag,
    bp_wet_dry_tag,
    bp_envfollow_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    constexpr uint  sr_taps_branch      = 32;
    constexpr uint  sr_taps_branch_frac = 16;
    constexpr float sr_cutoff           = 15000;
    constexpr float sr_kaiser_att_db    = 210;

    _resampler.reset (
      tgt_srate,
      pc.get_sample_rate(),
      sr_cutoff,
      sr_cutoff,
      sr_taps_branch,
      sr_taps_branch_frac,
      sr_kaiser_att_db,
      true,
      blocksize,
      6 * 1024);

    _n_spls_smoother.reset_coeffs (vec_set<4> (3.f), t_spl);
    _n_spls_smoother.reset_states();

    // get time info maximum buffer sizes and allocate
    _param.spls_x_beat = (60.f / pc.get_play_state().bpm) * (float) tgt_srate;
    initialize_buffer_related_parts();

    _filters.reset_states<lp_idx>();
    _filters.reset_states<hp_idx>();
    _filters.reset_states<bp_idx>();

    _tilt.reset_states();
    _transients.reset (tgt_srate);
    _transients.set (saike::transience::gainsmoothing_tag {}, 0.5f);

    using phase_type = decltype (_mod_lfo)::phase_type;
    phase_type::float_vec phases {0.f, 0.25f, 0.5f, 0.75f};
    _mod_lfo.reset();
    _mod_lfo.set_phase (phase_type {phase_type::normalized {}, phases});

    for (uint i = 0; i < _ap_lfo.size(); ++i) {
      _ap_lfo[i].reset();
      _ap_lfo[i].set_phase (phase_type {phase_type::normalized {}, phases});
      phases = vec_shuffle (phases, phases, 1, 2, 3, 0);
      phases += 0.01f;
      _ap_lfo[i].set_freq (vec_set<n_serial_diffusors> (0.247f), t_spl);
    }
#if DIFFUSE_DELAY_USE_SINC_TAPS
    _delay.reset_interpolator (0, false, 0.45f, 140.f);
#endif
#if DIFFUSE_DELAY_USE_THIRAN_DIFFUSORS
    for (uint i = 0; i < n_taps; ++i) {
      for (uint j = 0; j < n_serial_diffusors; ++j) {
        _diffusor[i][j].set_resync_delta (0.01f);
      }
    }
#endif

    // rms reset
    constexpr float rms_window_sec = 0.3f;
    _env.reset_coeffs<rms_dry_idx> (vec_set<4> (rms_window_sec), t_spl);
    _env.reset_coeffs<rms_wet_idx> (vec_set<4> (rms_window_sec), t_spl);
    _env.reset_coeffs<peakfollow_idx> (vec_set<4> (0.05f), t_spl);

    // hack to trigger intialization of some parameters
    _extpar.tilt_db           = 999.f;
    _extpar.bp_drive_db       = 999.f;
    _extpar.hp                = 999.f;
    _extpar.damp_ratio        = 999.f;
    _extpar.gain              = 999.f;
    _extpar.ducking_speed     = 999.f;
    _extpar.ducking_threshold = 999.f;

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    // avoid most of the initial smoother starting from time 0
    _n_spls_smoother.get_states()[0] = vec_set<n_taps> (_extpar.delay_spls);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    _resampler.process (outs, ins, samples, [=] (auto io) {
      process_block (io);
    });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  // GCD(44100,33600) = 2100. GCD(48000,33600) = 4800
  // for lower CPU: 31500: 6300 1500
  //                21000: 2100 3000
  static constexpr uint  tgt_srate          = 33600;
  static constexpr float t_spl              = 1.f / tgt_srate;
  static constexpr uint  blocksize          = 16;
  static constexpr uint  n_taps             = 4;
  static constexpr uint  n_serial_diffusors = 4;
  static constexpr uint  main_mod_samples   = 160; // bipolar (excursion 320)
  static constexpr uint  diffusor_mod_range = 48; // bipolar (excursion 96)
  using arith_type                          = float;
  using vec1_type                           = vec<arith_type, 1>;
  using vec_type                            = vec<arith_type, n_taps>;
  //----------------------------------------------------------------------------
  void initialize_buffer_related_parts()
  {
    // compute memory requirements
    // the headroom upwards is the minimum size of the delay (for)
    uint   headroom       = main_mod_samples + _delay.min_size_spls();
    double max_spls       = headroom + _param.spls_x_beat * max_t_beats;
    auto   fdn_size       = pow2_round_ceil ((uint) std::ceil (max_spls));
    _param.delay_spls_max = (double) (fdn_size);
    uint n_samples_delay  = _delay.n_required_elems (fdn_size, n_taps);

    auto                     allpass_sizes     = get_diffusor_delay_spls();
    uint                     n_samples_allpass = 0;
    uint const               n_allpasses       = allpass_sizes[0].size();
    std::array<uint, n_taps> row_sizes {};

    for (auto& row : allpass_sizes) {
      for (auto& delay : row) {
        headroom = diffusor_mod_range + _diffusor[0][0].min_size_spls();
        delay    = _diffusor[0][0].n_required_elems (
          pow2_round_ceil (delay + headroom), 1);
        n_samples_allpass += delay;
      }
    }
    // allocate
    _mem.clear();
    _mem.resize (n_samples_delay + n_samples_allpass);
    auto mem = xspan {_mem};

    // distribute
    _delay.reset (mem.cut_head (n_samples_delay), n_taps);
    for (uint i = 0; i < n_taps; ++i) {
      auto& line_sizes = allpass_sizes[i];
      for (uint j = 0; j < line_sizes.size(); ++j) {
        _diffusor[i][j].reset (mem.cut_head (line_sizes[j]), 1);
      }
    }
    assert (mem.size() == 0);
  }
  //----------------------------------------------------------------------------
  static std::array<std::array<uint, n_serial_diffusors>, n_taps>
  get_diffusor_delay_spls()
  {
    // some lines have the values of the Freeverb diffusor. Others come from
    // a schematic on the Gearslutz reverb subculture thread. Both work with
    // the allpasses at a gain of 0.5. Those were meant to be a starting point
    // but they already work due to the heavy modulation.
    return {
      {{{225u, 556u, 441u, 341u}},
       {{351u, 773u, 426u, 566u}},
       {{343u, 233u, 922u, 534u}},
       {{161u, 523u, 1171u, 1821u}}}};
  }
  //----------------------------------------------------------------------------
  double msec_to_spls (double msec)
  {
    return (double) tgt_srate * msec * 0.001;
  }
  //----------------------------------------------------------------------------
  void update_taps()
  {
    // main_mod samples wouldn't be necessary but this sets a minimum to have
    // the feedback time more or less stable.
    uint min_main_delay_spls
      = _delay.min_delay_spls() + blocksize + main_mod_samples;
    constexpr float max_delay_sec = 20.f;
    f32_x4          desync_spls_max {0.f, 161.803398f, 261.803f, 423.606f};

    auto delay_spls = std::max<float> (_extpar.delay_spls, min_main_delay_spls);
    f32_x4 spl_budget = vec_set<n_taps> (delay_spls - min_main_delay_spls);
    // desync
    auto desync_spls = spl_budget;
    spl_budget = vec_max (spl_budget - (desync_spls_max * _extpar.desync), 0.f);
    desync_spls -= spl_budget;
    // diffusors
    _param.diffusor_enable = 0;
    _param.diffusor_range  = sqrt (_extpar.diffusion);
    auto diffusor_n_spls   = get_diffusor_delay_spls();
    for (uint t = 0; t < n_taps; ++t) {
      for (uint d = 0; d < n_serial_diffusors; ++d) {
        float n_spls = diffusor_n_spls[t][d] * _param.diffusor_range;
        if (n_spls > spl_budget[t]) {
          continue;
        }
        spl_budget[t] -= n_spls;
        _param.diffusor_enable |= bit<u16> (t * n_taps + d);
      }
    }
    // main computation
    _param.delay_spls = spl_budget + (float) min_main_delay_spls;
    _param.fb_gain    = delay_get_feedback_gain_for_rt60_time (
      0.001f + _extpar.feedback * max_delay_sec,
      (float) tgt_srate,
      delay_spls - desync_spls);
    _param.max_hp_mod = log (1.f / _param.fb_gain[0]);
  }
  //----------------------------------------------------------------------------
  static constexpr std::array<float, 2> get_pan (
    float pan,
    float correction = 1.f)
  {
    std::array<float, 2> ret {};
    ret[0] = gcem::sin (M_PI_2 * (1.0f - pan)) * M_SQRT2 * correction;
    ret[1] = gcem::sin (M_PI_2 * pan) * M_SQRT2 * correction;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block (xspan<std::array<T, 2>> io)
  {
    auto const allpass_sizes = get_diffusor_delay_spls();

    std::array<vec_type, blocksize>               n_spls;
    array2d<arith_type, n_taps, blocksize>        tap_tail;
    array2d<arith_type, n_taps, blocksize>        tap_head;
    array2d<arith_type, 2, blocksize>             ducker_gain;
    array2d<float, n_serial_diffusors, blocksize> ap_spls;

    while (io.size()) {
      auto block = io.cut_head (std::min<uint> (io.size(), blocksize));
      // delay samples lfo
      auto depth      = (arith_type) (main_mod_samples * _extpar.mod_depth);
      auto delay_spls = _param.delay_spls;
      switch (_extpar.mod_mode) {
      // TODO enum instead of magic nums
      case 0:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i] = delay_spls;
          n_spls[i] += _mod_lfo.tick_filt_sample_and_hold() * depth;
        }
        break;
      case 1:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i] = delay_spls;
          n_spls[i] += _mod_lfo.tick_sine() * depth;
        }
        break;
      case 2:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i] = delay_spls;
          n_spls[i] += _mod_lfo.tick_triangle() * depth;
        }
        break;
      case 3:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i] = delay_spls;
          n_spls[i]
            += _mod_lfo.tick_trapezoid (vec_set<n_taps> (0.75f)) * depth;
        }
        break;
      default:
        assert (false);
        break;
      };
      // smoothing and clamping the delay in samples after modulation
      auto min_spls = vec_set<f32_x4> (_delay.min_delay_spls() + blocksize);
      for (uint i = 0; i < block.size(); ++i) {
        n_spls[i] = _n_spls_smoother.tick (n_spls[i]);
        n_spls[i] = vec_max (n_spls[i], min_spls);
      }
      // fill the tail samples, with feedback gain applied
      auto fb_gain = _param.fb_gain;
      for (uint t = 0; t < n_taps; ++t) {
        for (uint i = 0; i < block.size(); ++i) {
#if DIFFUSE_DELAY_USE_THIRAN_TAPS
          tap_tail[i][t] = _delay.get (n_spls[i][t], t, i)[0] * fb_gain[t];
#else
          tap_tail[i][t] = _delay.get (n_spls[i][t] - i, t)[0] * fb_gain[t];
#endif
        }
      }
      // tilt inputs
      for (uint i = 0; i < block.size(); ++i) {
        auto&  lr = block[i];
        f64_x2 ins {lr[0], lr[1]};
        auto   gain       = _ducker.tick (ins);
        ins               = _tilt.tick (ins);
        ducker_gain[i][0] = (arith_type) gain[0];
        ducker_gain[i][1] = (arith_type) gain[1];
        lr[0]             = (arith_type) ins[0];
        lr[1]             = (arith_type) ins[1];
      }

      // transient shaping
      _transients.process (block);
      // specific interleaving
      switch (_extpar.mode) {
      case m_stereo:
        stereo_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
      case m_ping_pong:
        pingpong_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
      case m_ping_pong_stereo:
        pingpong_stereo_interleaving<T> (
          tap_head.data(), block, tap_tail.data());
        break;
      case m_123:
      case m_132:
        x3_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
      case m_1223:
      case m_1234:
      case m_2413:
      case m_1324:
      case m_2314:
      case m_1423:
        x4_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
      case m_chorus:
        self_feed_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
      default:
        assert (false);
        break;
      }
      // diffusion by allpass
      auto gain            = make_vec (_extpar.diffusion * -0.5f);
      auto diffusor_enable = _param.diffusor_enable;
      auto diffusor_range  = _param.diffusor_range;

      for (uint t = 0; t < n_taps; ++t) {
        std::array<std::array<float, n_serial_diffusors>, blocksize> ap_spls;

        auto ap_spls_f = vec_cast<float> (vec_from_array (allpass_sizes[t]));
        for (uint i = 0; i < block.size(); ++i) {
          auto ap_lfo_f = _ap_lfo[t].tick_sine() * diffusor_mod_range;
          ap_lfo_f *= diffusor_range;
          ap_spls[i] = vec_to_array (ap_spls_f - ap_lfo_f);
        }

        for (uint d = 0; d < n_serial_diffusors; ++d) {
          bool enabled = !!(diffusor_enable & bit<u16> (t * n_taps + d));
          // As of now this is still run when disabled to ensure smooth
          // transitions, it might not be necessary.
          for (uint i = 0; i < block.size(); ++i) {
            auto v = allpass_fn::tick<vec1_type, float> (
              make_vec (tap_head[i][t]),
              ap_spls[i][d],
              gain,
              _diffusor[t][d])[0];
            tap_head[i][t] = enabled ? v : tap_head[i][t];
          }
        }
      }
      // diffusion by a 4-wide rotation matrix between taps
      for (uint i = 0; i < block.size(); ++i) {
        tap_head[i] = rotation_matrix<4>::tick<arith_type> (
          tap_head[i], _param.mtx_angle);
      }
      // Feedback FX
      float bp_drive = _extpar.bp_drive;
      float bp_wet   = _param.bp_wetdry;
      float bp_dry   = 1.f - abs (_param.bp_wetdry);
      float hp_gain  = fb_gain[0];
      hp_gain *= exp (_param.max_hp_mod * ((3.5f * _extpar.damp_ratio) - 2.5f));

      for (uint i = 0; i < block.size(); ++i) {
        auto taps = vec_from_array (tap_head[i]);
        // measuring feedback input power
        auto dry_rms = _env.tick<rms_dry_idx> (taps, envelope::rms_tag {});
        dry_rms      = vec_max (1e-30, dry_rms);

        // filter cutoff modulation by signal
        auto l     = block[i][0];
        auto r     = block[i][1];
        auto input = f32_x4 {(float) l, (float) r, (float) r, (float) r};
        auto input_env
          = _env.tick<peakfollow_idx> (input, envelope::rms_tag {});
        if ((bp_wet != 0.f) && ((_bp_update_spls & 15) == 0)) {
          auto freq = _param.bp_freqs;
          freq *= vec_exp (input_env * _extpar.bp_envfollow * 14.f);
          constexpr float filt_stability = 0.27f;
          freq
            = vec_min (freq, (float) (((tgt_srate / 2) - 1)) * filt_stability);
          _filters.reset_coeffs<bp_idx> (freq, get_scaled_reso (freq), t_spl);
        }
        ++_bp_update_spls;
        // Damp + HP/DC
        taps    = _filters.tick<hp_idx> (taps);
        auto lp = _filters.tick<lp_idx> (taps);
        auto hp = (taps - lp) * hp_gain; // hishelf
        taps    = lp + hp;

        // bping EQ FX
        if (bp_wet != 0.f) {
          auto wet = taps;
          wet *= bp_drive;
          wet = _filters.tick<bp_idx> (wet);
          wet *= bp_wet;
          taps *= bp_dry;
          taps += wet;
        }
        // measuring output power and gain riding the feedback gain
        auto wet_rms = _env.tick<rms_wet_idx> (taps, envelope::rms_tag {});
        wet_rms      = vec_max (1e-30, wet_rms);
        auto ratio   = dry_rms / (wet_rms);
        taps *= ratio;
        auto arr_x1 = vec1_array_wrap (vec_to_array (taps));
        _delay.push (arr_x1);
      }

      // specific output selection
      std::array<std::array<float, 2>, n_taps> tap_mul;

      constexpr auto pan_l    = get_pan (0.f);
      constexpr auto pan_l_x2 = get_pan (0.f, 2.f);
      constexpr auto pan_cl   = get_pan (0.333333f);
      constexpr auto pan_c    = get_pan (0.5f);
      constexpr auto pan_cr   = get_pan (0.666666f);
      constexpr auto pan_r    = get_pan (1.f);
      constexpr auto pan_r_x2 = get_pan (1.f, 2.f);
      constexpr auto zero     = get_pan (1.f, 0.f);
      switch (_extpar.mode) {
      case m_stereo:
      case m_ping_pong:
        tap_mul[0] = pan_l;
        tap_mul[1] = pan_l;
        tap_mul[2] = pan_r;
        tap_mul[3] = pan_r;
        break;
      case m_ping_pong_stereo:
        tap_mul[0] = pan_l_x2;
        tap_mul[1] = zero;
        tap_mul[2] = zero;
        tap_mul[3] = pan_r_x2;
        break;
      case m_123:
        tap_mul[0] = pan_l;
        tap_mul[1] = zero;
        tap_mul[2] = pan_c;
        tap_mul[3] = pan_r;
        break;
      case m_132:
        tap_mul[0] = pan_l;
        tap_mul[1] = zero;
        tap_mul[2] = pan_r;
        tap_mul[3] = pan_c;
        break;
      case m_1223:
        tap_mul[0] = pan_l;
        tap_mul[1] = pan_c;
        tap_mul[2] = pan_c;
        tap_mul[3] = pan_r;
        break;
      case m_1234:
        tap_mul[0] = pan_l;
        tap_mul[1] = pan_cl;
        tap_mul[2] = pan_cr;
        tap_mul[3] = pan_r;
        break;
      case m_2413:
        tap_mul[0] = pan_cl;
        tap_mul[1] = pan_r;
        tap_mul[2] = pan_l;
        tap_mul[3] = pan_cr;
        break;
      case m_1324:
        tap_mul[0] = pan_l;
        tap_mul[1] = pan_cr;
        tap_mul[2] = pan_cl;
        tap_mul[3] = pan_r;
        break;
      case m_2314:
        tap_mul[0] = pan_cl;
        tap_mul[1] = pan_cr;
        tap_mul[2] = pan_l;
        tap_mul[3] = pan_r;
        break;
      case m_1423:
        tap_mul[0] = pan_l;
        tap_mul[1] = pan_r;
        tap_mul[2] = pan_cl;
        tap_mul[3] = pan_cr;
        break;
      case m_chorus: {
        constexpr auto pan_l_2  = get_pan (0.f, 0.5f);
        constexpr auto pan_cl_2 = get_pan (0.333333f, 0.5f);
        constexpr auto pan_cr_2 = get_pan (0.666666f, 0.5f);
        constexpr auto pan_r_2  = get_pan (1.f, 0.5f);
        tap_mul[0]              = pan_l_2;
        tap_mul[1]              = pan_cl_2;
        tap_mul[2]              = pan_cr_2;
        tap_mul[3]              = pan_r_2;
      } break;
      default:
        assert (false);
        break;
      }
      // final tap acummulaton + gain
      auto main_gain = _param.main_gain;
      for (uint i = 0; i < block.size(); ++i) {
        // This is now tilted...
        block[i][0] = 0.f;
        block[i][1] = 0.f;
        for (uint t = 0; t < n_taps; ++t) {
          block[i][0] += tap_tail[i][t] * tap_mul[t][0];
          block[i][1] += tap_tail[i][t] * tap_mul[t][1];
        }
        block[i][0] *= main_gain * ducker_gain[i][0];
        block[i][1] *= main_gain * ducker_gain[i][1];
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void stereo_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    xspan<std::array<T, 2> const>         in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    // stereo with side ping pong
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][0];
      tap_head[i][1] = tap_tail[i][2];
      tap_head[i][2] = tap_tail[i][1];
      tap_head[i][3] = tap_tail[i][3];
      tap_head[i][0] += in[i][0];
      tap_head[i][3] += in[i][1];
      auto side = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][1] += side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void pingpong_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    xspan<std::array<T, 2> const>         in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][3];
      tap_head[i][1] = tap_tail[i][2];
      tap_head[i][2] = tap_tail[i][1];
      tap_head[i][3] = tap_tail[i][0];
      auto mid       = (arith_type) (in[i][0] + in[i][1]);
      auto side      = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][0] += mid;
      tap_head[i][1] += side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void pingpong_stereo_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    xspan<std::array<T, 2> const>         in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][1];
      tap_head[i][1] = tap_tail[i][0];
      tap_head[i][2] = tap_tail[i][3];
      tap_head[i][3] = tap_tail[i][2];
      tap_head[i][0] += in[i][0];
      tap_head[i][2] += in[i][1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void x3_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    xspan<std::array<T, 2> const>         in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][3];
      tap_head[i][1] = tap_tail[i][1];
      tap_head[i][2] = tap_tail[i][0];
      tap_head[i][3] = tap_tail[i][2];
      auto mid       = (arith_type) (in[i][0] + in[i][1]);
      tap_head[i][0] += mid;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void x4_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    xspan<std::array<T, 2> const>         in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][3];
      tap_head[i][1] = tap_tail[i][0];
      tap_head[i][2] = tap_tail[i][1];
      tap_head[i][3] = tap_tail[i][2];
      tap_head[i]
        = rotation_matrix<4>::tick<arith_type> (tap_head[i], _param.mtx_angle);
      auto mid  = (arith_type) (in[i][0] + in[i][1]);
      auto side = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][0] += mid + side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void self_feed_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    xspan<std::array<T, 2> const>         in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      auto mid = (arith_type) (in[i][0] + in[i][1]);
      for (uint j = 0; j < n_taps; ++j) {
        tap_head[i][j] = tap_tail[i][j];
      }
      for (uint j = 0; j < n_taps; ++j) {
        tap_head[i][j] += mid;
      }
    }
  }
  //----------------------------------------------------------------------------
  void update_damp()
  {
    constexpr float max_note  = 127.f; // 13289Hz
    constexpr float min_note  = 72.f;
    constexpr float bal_range = 20.f;

    float note = min_note + _extpar.damp_ratio * (max_note - min_note);
    _filters.reset_coeffs<lp_idx> (
      note_to_hzs (note, max_note, _extpar.freq_spread, bal_range), t_spl);
  }
  //----------------------------------------------------------------------------
  void update_bp()
  {
    constexpr float max_note  = constexpr_hz_to_midi_note (bp_max_hz);
    constexpr float bal_range = 6.f;
    constexpr float gain_mul  = 12.f;

    _param.bp_freqs
      = note_to_hzs (_extpar.bp_note, max_note, _extpar.freq_spread, bal_range);
    _param.bp_wetdry = _extpar.bp_wetdry * _extpar.bp_wetdry;
    _param.bp_wetdry *= _extpar.bp_wetdry >= 0.f ? 1.f : -1.f;
    // update the bp filter on the next block (it is done at mod rate)

    _bp_update_spls = 0;
  }
  //----------------------------------------------------------------------------
  f32_x4 get_scaled_reso (f32_x4 freq)
  {
    // high reso at high frequencies sounds ear pearcing, especially at this low
    // samplerate and low datatype resolution (for a nonlinear filter)
    constexpr float max_freq    = (float) ((tgt_srate - 1) / 2);
    constexpr float reso_att_hz = 1.f / max_freq;

    auto max_reso = vec_set<f32_x4> (_extpar.bp_reso);
    auto ratio    = freq * reso_att_hz;
    ratio         = vec_exp (-ratio * 2.f);
    return max_reso * ratio;
  }
  //----------------------------------------------------------------------------
  vec<arith_type, n_taps> note_to_hzs (
    float note,
    float max_note,
    float spread, //-1 to 1
    float spread_range_notes)
  {
    float diff    = abs (spread) * spread_range_notes * 0.5f;
    bool  reverse = (spread < 0.f);
    float current = note + (reverse ? diff : -diff);
    auto  step    = (diff * 2.f) / (n_taps - 1) * (reverse ? -1.f : 1.f);

    vec<float, n_taps> notes;
    for (uint i = 0; i < n_taps; ++i) {
      notes[i] = current;
      current += step;
    }
    return midi_note_to_hz (vec_min (notes, max_note));
  }
  //----------------------------------------------------------------------------
  struct external_parameters {
    uint  mode;
    uint  mod_mode;
    float delay_spls;
    float feedback;
    float diffusion;
    float gain;
    float ducking_threshold;
    float ducking_speed;
    float mod_freq;
    float mod_depth;
    float mod_spread;
    float damp_ratio;
    float freq_spread;
    float tilt_db;
    float desync;
    float transients;
    float hp;
    float bp_freq;
    float bp_envfollow;
    float bp_drive_db;
    float bp_drive;
    float bp_wetdry;
    float bp_note;
    float bp_reso;
  };
  //----------------------------------------------------------------------------
  struct internal_parameters {
    f32_x4               bp_freqs;
    f32_x4               delay_spls;
    f32_x4               fb_gain;
    array2d<float, 2, 1> mtx_angle;
    double               delay_spls_max;
    double               spls_x_beat;
    float                bp_wetdry;
    float                diffusor_range;
    float                main_gain;
    float                max_hp_mod;
    u16                  diffusor_enable;
  };
//----------------------------------------------------------------------------
#if DIFFUSE_DELAY_USE_THIRAN_DIFFUSORS
  std::array<
    std::array<modulable_thiran1_delay_line<vec1_type>, n_serial_diffusors>,
    n_taps>
    _diffusor;
#else
  std::array<
    std::array<
      interpolated_delay_line<vec1_type, catmull_rom_interp>,
      n_serial_diffusors>,
    n_taps>
    _diffusor;
#endif
  //----------------------------------------------------------------------------
  // in case sinc interpolation is enabled, it is nice to have the tables
  // aligned to cache line boundaries.
  template <class T>
  using mem_vector = std::vector<T, overaligned_allocator<T, 128>>;

  uint                           _bp_update_spls {};
  double                         _gr_prev {};
  external_parameters            _extpar {};
  internal_parameters            _param {};
  block_resampler<arith_type, 2> _resampler {};
#if DIFFUSE_DELAY_USE_THIRAN_TAPS
  modulable_thiran1_delay_line<vec1_type, 4> _delay {};
#elif DIFFUSE_DELAY_USE_SINC_TAPS
  interpolated_delay_line<vec1_type, sinc_interp<16, 128>> _delay {};
#else
  interpolated_delay_line<vec1_type, catmull_rom_interp> _delay {};
#endif
  mem_vector<vec1_type>                       _mem {};
  part_class_array<onepole_smoother, f32_x4>  _n_spls_smoother {};
  part_class_array<tilt_eq, f64_x2>           _tilt {};
  saike::transience                           _transients;
  lfo<n_taps>                                 _mod_lfo;
  std::array<lfo<n_serial_diffusors>, n_taps> _ap_lfo;
  ducker<f64_x2>                              _ducker;
  enum { rms_dry_idx, rms_wet_idx, peakfollow_idx };
  part_classes<mp_list<envelope, envelope, envelope>, f32_x4, false> _env;
  enum { lp_idx, hp_idx, bp_idx };
  part_classes<
    mp_list<onepole_lowpass, mystran_dc_blocker, saike::ms20_bandpass>,
    vec<arith_type, n_taps>,
    false>
    _filters;
};
//------------------------------------------------------------------------------
} // namespace artv
