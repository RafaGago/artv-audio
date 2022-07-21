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
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/interpolation/stateless.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/third_party/saike/transience.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#define DIFFUSE_DELAY_USE_THIRAN_TAPS 1
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
    delay_line_updated();
  }

  static constexpr auto get_parameter (sixteenths_tag)
  {
    return float_param ("sixteenths", 0.1f, max_t_beats * 16, 6.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    v = v * 0.01f;
    if (v == _extpar.feedback) {
      return;
    }
    _extpar.feedback = v;
    delay_line_updated();
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("%", 0.f, 100.f, 25.f, 0.001f);
  }
  //----------------------------------------------------------------------------
  struct diffusion_tag {};
  void set (diffusion_tag, float v)
  {
    v *= 0.01f;
    // v *= v;
    v *= -0.5f; // cut at 0.5f for the allpass
    if (v == _extpar.diffusion) {
      return;
    }
    _extpar.diffusion = v;
  }

  static constexpr auto get_parameter (diffusion_tag)
  {
    return float_param ("%", 0.f, 100.f, 20.f, 0.01f);
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
    _ducker.set_speed (vec_set<double_x2> (v * v), tgt_srate);
  }

  static constexpr auto get_parameter (ducking_speed_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
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
    return float_param ("dB", -50.f, 12.f, 0.f, 0.2f);
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
    _mod_lfo.set_freq (vec_set<n_taps> (v), (float) tgt_srate);
  }

  static constexpr auto get_parameter (mod_freq_tag)
  {
    return float_param ("Hz", 0.f, 24.f, 0.17f, 0.001f, 0.35f);
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
  struct damp_freq_tag {};
  void set (damp_freq_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.damp_note_ratio) {
      return;
    }
    _extpar.damp_note_ratio = v;
    update_damp();
  }

  static constexpr auto get_parameter (damp_freq_tag)
  {
    return float_param ("%", 0.f, 100.f, 75.f, 0.01f);
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
    update_peak();
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
      (double) tgt_srate);
  }

  static constexpr auto get_parameter (tilt_db_tag)
  {
    return float_param ("dB", -16.f, 16.f, 0.f, 0.1f);
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
    _filters.reset_coeffs<hp_idx> (
      vec_set<n_taps> (4.f + v * 396.f), (float) tgt_srate);
  }

  static constexpr auto get_parameter (hipass_tag)
  {
    return float_param ("%", 0.f, 100.f, 25.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct peak_drive_tag {};
  void set (peak_drive_tag, float v)
  {
    // -30dB to +30dB range
    constexpr float dbrange = 30.f;
    v *= (0.01 * dbrange * 2.f);
    v -= dbrange;
    if (v == _extpar.peak_drive_db) {
      return;
    }
    _extpar.peak_drive_db = v;
    _extpar.peak_drive    = db_to_gain (v);
  }

  static constexpr auto get_parameter (peak_drive_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct peak_freq_tag {};
  void set (peak_freq_tag, float v)
  {
    if (v == _extpar.peak_note) {
      return;
    }
    _extpar.peak_note = v;
    update_peak();
  }

  static constexpr float peak_min_hz = 100.;
  static constexpr float peak_max_hz = 5000.;

  static constexpr auto get_parameter (peak_freq_tag)
  {
    return frequency_parameter (peak_min_hz, peak_max_hz, 440.0);
  }
  //----------------------------------------------------------------------------
  struct peak_gain_tag {};
  void set (peak_gain_tag, float v)
  {
    v *= 0.01;
    if (v == _extpar.peak_gain) {
      return;
    }
    _extpar.peak_gain = v;
    update_peak();
  }

  static constexpr auto get_parameter (peak_gain_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    feedback_tag,
    gain_tag,
    sixteenths_tag,
    damp_freq_tag,
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
    peak_drive_tag,
    peak_freq_tag,
    peak_gain_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    constexpr uint  sr_taps_branch      = 32;
    constexpr uint  sr_taps_branch_frac = 16;
    constexpr float sr_cutoff           = 15000;
    constexpr float sr_kaiser_att_db    = 210;

    _resampler.reset (
      pc.get_sample_rate(),
      tgt_srate,
      sr_cutoff,
      sr_cutoff,
      sr_taps_branch,
      sr_taps_branch_frac,
      sr_kaiser_att_db,
      true,
      blocksize,
      6 * 1024);

    _n_spls_smoother.reset_coeffs (vec_set<4> (1.f), tgt_srate);
    _n_spls_smoother.reset_states();

    // get time info maximum buffer sizes and allocate
    _param.spls_x_beat = pc.get_samples_per_beat();
    initialize_buffer_related_parts();

    _filters.reset_states<lp_idx>();
    _filters.reset_states<hp_idx>();
    _filters.reset_states<peak_idx>();

    _tilt.reset_states();
    _transients.reset (tgt_srate);
    _transients.set (saike::transience::gainsmoothing_tag {}, 0.5f);

    using phase_type = decltype (_mod_lfo)::phase_type;
    using value_type = decltype (_mod_lfo)::value_type;

    value_type phases {0.f, 0.25f, 0.5f, 0.75f};

    _mod_lfo.reset();
    _mod_lfo.set_phase (phase_type {phases, phase_type::normalized {}});

    for (uint i = 0; i < _ap_lfo.size(); ++i) {
      _ap_lfo[i].reset();
      _ap_lfo[i].set_phase (phase_type {phases, phase_type::normalized {}});
      phases = vec_shuffle (phases, phases, 1, 2, 3, 0);
      phases += 0.01f;
      _ap_lfo[i].set_freq (
        vec_set<n_serial_diffusors> (0.247f), (float) tgt_srate);
    }
#if DIFFUSE_DELAY_USE_THIRAN_TAPS
    _delay.set_resync_delta (0.0);
#endif
#if DIFFUSE_DELAY_USE_THIRAN_DIFFUSORS
    for (uint i = 0; i < n_taps; ++i) {
      for (uint j = 0; j < n_serial_diffusors; ++j) {
        _diffusor[i][j].set_resync_delta (0.01f);
      }
    }
#endif

    // rms reset
    constexpr float rms_window_sec = 0.3;
    _rms.reset_coeffs<rms_dry_idx> (vec_set<4> (rms_window_sec), tgt_srate);
    _rms.reset_coeffs<rms_wet_idx> (vec_set<4> (rms_window_sec), tgt_srate);

    // hack to trigger intialization of some parameters
    _extpar.tilt_db           = 999.f;
    _extpar.peak_drive_db     = 999.f;
    _extpar.hp                = 999.f;
    _extpar.damp_note_ratio   = 999.f;
    _extpar.gain              = 999.f;
    _extpar.ducking_speed     = 999.f;
    _extpar.ducking_threshold = 999.f;

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
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
  static constexpr uint tgt_srate          = 33600;
  static constexpr uint blocksize          = 16;
  static constexpr uint n_taps             = 4;
  static constexpr uint n_serial_diffusors = 4;
  static constexpr uint max_mod_samples    = 192; // unipolar
  static constexpr uint diffusor_mod_range = 48; // bipolar
  using arith_type                         = float;
  using vec1_type                          = vec<arith_type, 1>;
  using vec_type                           = vec<arith_type, n_taps>;
  //----------------------------------------------------------------------------
  void initialize_buffer_related_parts()
  {
    // compute memory requirements
    double max_spls = max_mod_samples + diffusor_mod_range;
    max_spls += _param.spls_x_beat * max_t_beats;
    auto fdn_size         = pow2_round_ceil ((uint) std::ceil (max_spls));
    _param.delay_spls_max = (double) (fdn_size);
    uint n_samples_delay  = _delay.n_required_elems (fdn_size, n_taps);

    auto                     allpass_sizes     = get_diffusor_delay_spls();
    uint                     n_samples_allpass = 0;
    uint const               n_allpasses       = allpass_sizes[0].size();
    std::array<uint, n_taps> row_sizes {};

    for (auto& row : allpass_sizes) {
      for (auto& delay : row) {
        delay = _diffusor[0][0].n_required_elems (
          pow2_round_ceil (delay + diffusor_mod_range), 1);
        n_samples_allpass += delay;
      }
    }
    // allocate
    _mem.clear();
    _mem.resize (n_samples_delay + n_samples_allpass);
    auto mem = make_crange (_mem);

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
    // but I already liked them from the start.
    return {
      {{{225u, 556u, 441u, 341u}},
       {{351u, 773u, 426u, 566u}},
       {{343u, 233u, 922u, 534u}},
       {{161u, 523u, 1171u, 1821u}}}};
  }
  //----------------------------------------------------------------------------
  static std::array<uint, n_taps> get_diffusor_delay_total_spls()
  {
    auto                     spls = get_diffusor_delay_spls();
    std::array<uint, n_taps> ret;
    for (uint i = 0; i < ret.size(); ++i) {
      uint sum = 0;
      for (auto v : spls[i]) {
        sum += v;
      }
      ret[i] = sum;
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  double msec_to_spls (double msec)
  {
    return (double) tgt_srate * msec * 0.001;
  }
  //----------------------------------------------------------------------------
  void delay_line_updated()
  {
    _param.fb_gain = delay_get_feedback_gain_for_time (
      _extpar.feedback * 20.f,
      -60.f,
      (float) tgt_srate,
      vec_set<1> (_extpar.delay_spls > 0.f ? _extpar.delay_spls : 0.1f))[0];
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
  void process_block (crange<std::array<T, 2>> io)
  {
    auto const diffusor_correction
      = vec_cast<arith_type> (vec_from_array (get_diffusor_delay_total_spls()));

    auto const allpass_sizes = get_diffusor_delay_spls();

    std::array<vec_type, blocksize>                              n_spls;
    std::array<std::array<arith_type, n_taps>, blocksize>        tap_tail;
    std::array<std::array<arith_type, n_taps>, blocksize>        tap_head;
    std::array<std::array<float, n_serial_diffusors>, blocksize> ap_spls;

    while (io.size()) {
      auto block = io.cut_head (std::min<uint> (io.size(), blocksize));

      // delay samples smoothed
      auto n_spls_readonce = vec_set<n_taps> (_extpar.delay_spls);
      auto desync          = _extpar.desync;
      for (uint i = 0; i < block.size(); ++i) {

        n_spls[i] = n_spls_readonce;
        n_spls[i] -= diffusor_correction;
        float_x4 const desync_spls_max {0.f, 161.803398f, 261.803f, 423.606f};
        n_spls[i] -= desync_spls_max * desync;
      }
      // delay samples lfo
      auto depth = (arith_type) (max_mod_samples * _extpar.mod_depth);
      switch (_extpar.mod_mode) {
      // TODO enum instead of magic nums
      case 0:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i] += _mod_lfo.tick_filt_sample_and_hold() * depth;
        }
        break;
      case 1:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i] += _mod_lfo.tick_sine() * depth;
        }
        break;
      case 2:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i] += _mod_lfo.tick_triangle() * depth;
        }
        break;
      case 3:
        for (uint i = 0; i < block.size(); ++i) {
          n_spls[i]
            += _mod_lfo.tick_trapezoid (vec_set<n_taps> (0.75f)) * depth;
        }
        break;
      default:
        assert (false);
        break;
      };
      // smoothing and clamping the delay in samples after modulation
      for (uint i = 0; i < block.size(); ++i) {
        n_spls[i] = _n_spls_smoother.tick (n_spls[i]);
#if DIFFUSE_DELAY_USE_THIRAN_TAPS
        // 2 which is the number of states of a second order filter (Thiran2)
        float const interp_headroom_spls
          = _delay.n_interp_states + _delay.n_warmup_spls;
#else
        float const interp_headroom_spls = catmull_rom_interp::n_points;
#endif
        auto min_spls = vec_set<n_taps> (interp_headroom_spls + blocksize);
        n_spls[i]     = n_spls[i] >= min_spls ? n_spls[i] : min_spls;
        n_spls[i] -= vec_set<n_taps> ((arith_type) i);
      }
      // fill the tail samples, with feedback gain applied
      auto fb_gain = _param.fb_gain;
      for (uint t = 0; t < n_taps; ++t) {
        for (uint i = 0; i < block.size(); ++i) {
          tap_tail[i][t] = _delay.get (n_spls[i][t], t)[0] * fb_gain;
        }
      }
      // transient shaping
      _transients.process (block);
      // tilt inputs
      for (uint i = 0; i < block.size(); ++i) {
        auto&     lr = block[i];
        double_x2 tilted {lr[0], lr[1]};
        tilted = _tilt.tick (tilted);
        lr[0]  = (arith_type) tilted[0];
        lr[1]  = (arith_type) tilted[1];
      }
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
      // diffusion
      auto gain = make_vec (_extpar.diffusion);
      for (uint t = 0; t < n_taps; ++t) {
        std::array<std::array<float, n_serial_diffusors>, blocksize> ap_spls;

        auto ap_spls_f = vec_cast<float> (vec_from_array (allpass_sizes[t]));
        for (uint i = 0; i < block.size(); ++i) {
          auto ap_lfo_f = _ap_lfo[t].tick_sine() * diffusor_mod_range;
          ap_spls[i]    = vec_to_array (ap_spls_f - ap_lfo_f);
        }

        for (uint d = 0; d < n_serial_diffusors; ++d) {
          for (uint i = 0; i < block.size(); ++i) {
            tap_head[i][t] = allpass_fn::tick<vec1_type, float> (
              make_vec (tap_head[i][t]),
              ap_spls[i][d],
              gain,
              _diffusor[t][d])[0];
          }
        }
      }
      // Feedback FX
      float peak_drive     = _extpar.peak_drive;
      float peak_drive_inv = 1. / peak_drive;
      for (uint i = 0; i < block.size(); ++i) {
        auto taps = vec_from_array (tap_head[i]);
        // measuring input power
        auto dry_rms = _rms.tick<rms_dry_idx> (taps, envelope::rms_tag {});
        dry_rms      = vec_max (1e-30, dry_rms);
        // Damp + HP/DC
        taps = _filters.tick<lp_idx> (taps);
        taps = _filters.tick<hp_idx> (taps);
        // Peaking EQ FX
        auto wet = _filters.tick<peak_idx> (taps);
        wet *= peak_drive;
        wet = wet / vec_sqrt (1.f + wet * wet);
        wet *= peak_drive_inv;
        taps += wet;
        // measuring output power and gain riding the feedback gain
        auto wet_rms = _rms.tick<rms_wet_idx> (taps, envelope::rms_tag {});
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
        auto ducking_gain = _ducker.tick (double_x2 {block[i][0], block[i][1]});
        block[i][0]       = 0.f;
        block[i][1]       = 0.f;
        for (uint t = 0; t < n_taps; ++t) {
          block[i][0] += tap_tail[i][t] * tap_mul[t][0];
          block[i][1] += tap_tail[i][t] * tap_mul[t][1];
        }
        auto out = double_x2 {block[i][0], block[i][1]};
        out *= main_gain * ducking_gain;
        block[i][0] = (T) out[0];
        block[i][1] = (T) out[1];
      }
    }
  }
  //----------------------------------------------------------------------------
  auto get_fdn4_angle()
  {
    // _extpar.diffusion goes from 0 to 0.5
    float weight = _extpar.diffusion * (2. * (M_PI / 4));
    std::array<std::array<float, 2>, 1> angle;
    angle[0][0] = cos (weight);
    angle[0][1] = sqrt (1.f - angle[0][0] * angle[0][0]);
    return angle;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void stereo_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    crange<std::array<T, 2> const>        in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    // stereo with side ping pong
    auto angle = get_fdn4_angle();
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][0];
      tap_head[i][1] = tap_tail[i][2];
      tap_head[i][2] = tap_tail[i][1];
      tap_head[i][3] = tap_tail[i][3];
      // Extra FDN diffusion.
      tap_head[i] = rotation_matrix<4>::tick<arith_type> (tap_head[i], angle);
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
    crange<std::array<T, 2> const>        in, // inputs
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
    crange<std::array<T, 2> const>        in, // inputs
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
    crange<std::array<T, 2> const>        in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][3];
      tap_head[i][1] = 0.f; // disabled
      tap_head[i][2] = tap_tail[i][0] + tap_tail[i][1];
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
    crange<std::array<T, 2> const>        in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][3];
      tap_head[i][1] = tap_tail[i][0];
      tap_head[i][2] = tap_tail[i][1];
      tap_head[i][3] = tap_tail[i][2];
      auto mid       = (arith_type) (in[i][0] + in[i][1]);
      auto side      = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][0] += mid + side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void self_feed_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be
                                              // inserted
    crange<std::array<T, 2> const>        in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    auto angle = get_fdn4_angle();
    for (uint i = 0; i < in.size(); ++i) {
      auto mid = (arith_type) (in[i][0] + in[i][1]);
      for (uint j = 0; j < n_taps; ++j) {
        tap_head[i][j] = tap_tail[i][j];
      }
      tap_head[i] = rotation_matrix<4>::tick<arith_type> (tap_head[i], angle);
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

    float note = min_note + _extpar.damp_note_ratio * (max_note - min_note);
    _filters.reset_coeffs<lp_idx> (
      note_to_hzs (note, max_note, _extpar.freq_spread, bal_range),
      (float) tgt_srate);
  }
  //----------------------------------------------------------------------------
  void update_peak()
  {
    constexpr float min_note    = constexpr_hz_to_midi_note (peak_min_hz);
    constexpr float gr_note     = constexpr_hz_to_midi_note (2000.);
    constexpr float max_note    = constexpr_hz_to_midi_note (peak_max_hz);
    constexpr float note_weight = 1.f / (max_note - gr_note);
    constexpr float bal_range   = 4.f;

    // don't reduce the gain at higher frequencies, it is disgusting
    auto  note    = _extpar.peak_note;
    auto  absgain = abs (_extpar.peak_gain);
    float peak_reduction
      = (note < gr_note) ? 0.f : (note - gr_note) * note_weight;
    float max_db = (7.f - peak_reduction * absgain * 2.f);

    _filters.reset_coeffs<peak_idx> (
      note_to_hzs (note, max_note, _extpar.freq_spread, bal_range),
      vec_set<4> (0.2f + absgain * 0.6f), // Q
      vec_set<4> (_extpar.peak_gain * max_db), // dB
      (float) tgt_srate,
      bell_bandpass_tag {});
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
    float damp_note_ratio;
    float freq_spread;
    float tilt_db;
    float desync;
    float transients;
    float hp;
    float peak_freq;
    float peak_drive_db;
    float peak_drive;
    float peak_gain;
    float peak_note;
  };
  //----------------------------------------------------------------------------
  struct internal_parameters {
    double delay_spls_max;
    double spls_x_beat;
    float  fb_gain;
    float  main_gain;
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
  double                         _gr_prev {};
  external_parameters            _extpar {};
  internal_parameters            _param {};
  block_resampler<arith_type, 2> _resampler {};
#if DIFFUSE_DELAY_USE_THIRAN_TAPS
  modulable_thiran1_delay_line<vec1_type> _delay {};
#else
  interpolated_delay_line<vec1_type, catmull_rom_interp> _delay {};
#endif
  std::vector<vec1_type>                       _mem {};
  part_class_array<onepole_smoother, float_x4> _n_spls_smoother {};
  part_class_array<tilt_eq, double_x2>         _tilt {};
  saike::transience                            _transients;
  lfo<n_taps>                                  _mod_lfo;
  std::array<lfo<n_serial_diffusors>, n_taps>  _ap_lfo;
  ducker<double_x2>                            _ducker;
  enum { rms_dry_idx, rms_wet_idx };
  part_classes<mp_list<envelope, envelope>, float_x4, false> _rms;
  enum { lp_idx, hp_idx, peak_idx };
  part_classes<
    mp_list<onepole_lowpass, mystran_dc_blocker, andy::svf>,
    vec<arith_type, n_taps>,
    false>
    _filters;
};
//------------------------------------------------------------------------------
} // namespace artv
