#pragma once

#include <gcem.hpp>
#include <numeric>
#include <vector>

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/diffusion_matrix.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/reverb_tools.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/interpolation/stateless.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

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
    return float_param ("sixteenths", 0.1f, max_t_beats * 16, 12.f, 0.01f);
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
    return float_param ("%", 0.f, 100.f, 25.f, 0.01f);
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
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct ducking_threshold_tag {};
  void set (ducking_threshold_tag, float v) { _extpar.ducking_threshold = v; }

  static constexpr auto get_parameter (ducking_threshold_tag)
  {
    return float_param ("dB", -90.f, 20.f, 20.f, 0.1f, 1.5f);
  }
  //----------------------------------------------------------------------------
  struct ducking_speed_tag {};
  void set (ducking_speed_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.ducking_speed) {
      return;
    }
    _extpar.ducking_speed = v;
    _ducker_follow.reset_coeffs (
      vec_set<1> (0.03f),
      vec_set<1> (0.01f + 0.1f * v),
      (float) sr_target_freq);
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
    return float_param ("dB", -50.f, 12.f, -20.f, 0.2f);
  }
  //----------------------------------------------------------------------------
  struct mode_param_tag {};
  void set (mode_param_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.mode_param) {
      return;
    }
    _extpar.mode_param = v;
  }

  static constexpr auto get_parameter (mode_param_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
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
    m_stereo2,
    m_ping_pong,
    m_123,
    m_132,
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
        "Stereo 2",
        "Ping-Pong",
        "P1-P2-P3 (L-C-R)",
        "P1-P3-P2 (L-R-C)",
        "P1-P2-P3-P4",
        "P2-P4-P1-P3",
        "P1-P3-P2-P4",
        "P2-P3-P1-P4",
        "P1-P4-P2-P3",
        "Chorus-friendy"),
  }

  //----------------------------------------------------------------------------
  struct mod_freq_tag {};
  void set (mod_freq_tag, float v)
  {
    if (v == _extpar.mod_freq) {
      return;
    }
    _extpar.mod_freq = v;
    _mod_lfo.set_freq (vec_set<n_taps> (v), (float) sr_target_freq);
  }

  static constexpr auto get_parameter (mod_freq_tag)
  {
    return float_param ("Hz", 0.f, 4.f, 0.17f, 0.01f, 0.5f);
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
    return float_param ("%", 0.f, 100.f, 25.f, 0.01f);
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
    if (v == _extpar.damp_freq) {
      return;
    }
    constexpr float max_note   = 127; // 13289Hz
    constexpr float min_note   = 60; // 261Hz
    constexpr float note_range = max_note - min_note;

    _extpar.damp_freq = v;
    float freq        = midi_note_to_hz (min_note + v * note_range);
    _filters.reset_coeffs<lp_idx> (
      vec_set<4> ((float) freq), (float) sr_target_freq);
  }

  static constexpr auto get_parameter (damp_freq_tag)
  {
    return float_param ("%", 0.f, 100.f, 50.f, 0.01f);
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
      (double) sr_target_freq);
  }

  static constexpr auto get_parameter (tilt_db_tag)
  {
    return float_param ("dB", -16.f, 16.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct stereo_tag {};
  void set (stereo_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.stereo) {
      return;
    }
    _extpar.stereo = v;
  }

  static constexpr auto get_parameter (stereo_tag)
  {
    return float_param ("%", 0.f, 100.f, 100.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    feedback_tag,
    gain_tag,
    sixteenths_tag,
    damp_freq_tag,
    tilt_db_tag,
    mode_tag,
    mode_param_tag,
    diffusion_tag,
    stereo_tag,
    mod_mode_tag,
    mod_freq_tag,
    mod_depth_tag,
    ducking_speed_tag,
    ducking_threshold_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    // GCD(44100,33600) = 2100. GCD(48000,33600) = 4800
    constexpr uint  sr_taps_branch      = 32;
    constexpr uint  sr_taps_branch_frac = 16;
    constexpr float sr_cutoff           = 12000;
    constexpr float sr_kaiser_att_db    = 210;

    _resampler.reset (
      pc.get_sample_rate(),
      sr_target_freq,
      sr_cutoff,
      sr_cutoff,
      sr_taps_branch,
      sr_taps_branch_frac,
      sr_kaiser_att_db,
      true,
      blocksize,
      6 * 1024);

    _n_spls_smoother.reset_coeffs (vec_set<1> (1.f), sr_target_freq);
    _n_spls_smoother.reset_states();

    // get time info maximum buffer sizes and allocate
    _param.spls_x_beat = pc.get_samples_per_beat();
    initialize_buffer_related_parts();

    _filters.reset_coeffs<dc_idx> (
      vec_set<n_taps> (1.f), (float) sr_target_freq);

    _filters.reset_states<lp_idx>();
    _filters.reset_states<dc_idx>();

    _tilt.reset_states();

    using phase_type = decltype (_mod_lfo)::phase_type;
    using value_type = decltype (_mod_lfo)::value_type;

    _mod_lfo.reset();
    _mod_lfo.set_phase (phase_type {
      value_type {0.f, 0.25f, 0.5f, 0.75f}, phase_type::normalized {}});

    _ap_lfo.reset();
    _ap_lfo.set_phase (phase_type {
      value_type {0.f, 0.25f, 0.5f, 0.75f}, phase_type::normalized {}});
    _ap_lfo.set_freq (vec_set<n_taps> (0.1f), (float) sr_target_freq);

    // hack to trigger intialization of the tilt filter before so a change is
    // detected on "mp_for_each"
    _extpar.tilt_db = 1.f;

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
  static constexpr uint blocksize       = 16;
  static constexpr uint n_taps          = 4;
  static constexpr uint max_mod_samples = 500;
  using arith_type                      = float;
  using vec1_type                       = vec<arith_type, 1>;
  using vec_type                        = vec<arith_type, n_taps>;
  //----------------------------------------------------------------------------
  void initialize_buffer_related_parts()
  {
    // compute memory requirements
    double max_spls = max_mod_samples;
    max_spls += _param.spls_x_beat * max_t_beats;
    auto fdn_size         = pow2_round_ceil ((uint) std::ceil (max_spls));
    _param.delay_spls_max = (double) (fdn_size);
    uint n_samples_delay  = _delay.n_required_elems (n_taps, fdn_size);

    auto allpass_sizes     = get_diffusor_delay_spls();
    uint n_samples_allpass = 0;
    for (auto& row : allpass_sizes) {
      for (auto& delay : row) {
        delay = pow2_round_ceil (delay);
        n_samples_allpass += delay;
      }
    }

    // allocate
    _mem.clear();
    _mem.resize (n_samples_delay + n_samples_allpass);
    auto mem = make_crange (_mem);

    // distribute
    _delay.reset (mem.cut_head (n_samples_delay), n_taps);

    for (uint i = 0; i < allpass_sizes.size(); ++i) {
      auto& line_sizes = allpass_sizes[i];
      for (uint j = 0; j < line_sizes.size(); ++j) {
        _diffusor[i][j].reset (mem.cut_head (allpass_sizes[i][j]));
      }
    }
  }
  //----------------------------------------------------------------------------
  static std::array<std::array<uint, 4>, n_taps> get_diffusor_delay_spls()
  {
    // some lines have the values of the Freeverb diffusor. Others come from
    // a schematic on the Gearslutz reverb subculture thread. Both work with
    // the allpasses at a gain of 0.5. Those were meant to be a starting point
    // but I already liked them from the start.
    return {
      {{{225u, 556u, 441u, 341u}},
       {{161, 523u, 1171u, 1821u}},
       {{225u, 556u, 441u, 341u}},
       {{161, 523u, 1171u, 1821u}}}};
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
  static constexpr uint sr_target_freq = 33600;
  //----------------------------------------------------------------------------
  double msec_to_spls (double msec)
  {
    return (double) sr_target_freq * msec * 0.001;
  }
  //----------------------------------------------------------------------------
  void delay_line_updated()
  {
    _param.fb_gain = delay_get_feedback_gain_for_time (
      _extpar.feedback * 20.f,
      -60.f,
      (float) sr_target_freq,
      vec_set<1> (_extpar.delay_spls > 0.f ? _extpar.delay_spls : 0.1f))[0];
  }
  //----------------------------------------------------------------------------
  float get_ducker_gain (float in)
  {
    // https://www.musicdsp.org/en/latest/Effects/204-simple-compressor-class-c.html

    constexpr float ratio = 9.f;

    in          = gain_to_db (fabs (in), -240.f); // convert linear -> dB
    float delta = in - _extpar.ducking_threshold;
    delta       = std::max (delta, 0.f);
    float env   = _ducker_follow.tick (vec_set<1> (delta))[0];
    float gr    = env * (ratio - 1.f);
    return db_to_gain (-gr);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block (crange<std::array<T, 2>> io)
  {
    auto const diffusor_correction
      = vec_cast<arith_type> (vec_from_array (get_diffusor_delay_total_spls()));
    auto const allpass_sizes = get_diffusor_delay_spls();

    while (io.size()) {
      auto block = io.cut_head (std::min<uint> (io.size(), blocksize));

      std::array<vec_type, blocksize>                       n_spls;
      std::array<std::array<arith_type, n_taps>, blocksize> tap_head;
      std::array<std::array<arith_type, n_taps>, blocksize> tap_tail;

      // delay samples smoothed
      auto del_spls = vec_set<1> (_extpar.delay_spls);
      for (uint i = 0; i < block.size(); ++i) {
        auto smooth = _n_spls_smoother.tick (del_spls)[0];
        n_spls[i]   = vec_set<n_taps> ((arith_type) smooth);
        n_spls[i] -= diffusor_correction;
      }
      // delay samples lfo
      auto mode  = _extpar.mod_mode;
      auto depth = (arith_type) (max_mod_samples * _extpar.mod_depth);
      for (uint i = 0; i < block.size(); ++i) {
        switch (mode) {
        // TODO enum instead of magic nums
        case 0:
          n_spls[i] += _mod_lfo.tick_filt_sample_and_hold() * depth;
          break;
        case 1:
          n_spls[i] += _mod_lfo.tick_sine() * depth;
          break;
        case 2:
          n_spls[i] += _mod_lfo.tick_triangle() * depth;
          break;
        case 3:
          n_spls[i]
            += _mod_lfo.tick_trapezoid (vec_set<n_taps> (0.75f)) * depth;
          break;
        default:
          assert (false);
          break;
        };
      }
      // clamping the delay in samples after modulation
      for (uint i = 0; i < block.size(); ++i) {
        // 2 which is the number of states of a second order filter (Thiran2)
        auto min_spls = vec_set<n_taps> (2.f + blocksize);
        n_spls[i]     = n_spls[i] >= min_spls ? n_spls[i] : min_spls;
        n_spls[i] -= vec_set<n_taps> ((arith_type) i);
      }
      // fill the tail samples, with feedback gain applied
      for (uint i = 0; i < block.size(); ++i) {
        std::array<vec1_type, n_taps> tailv;
        auto                          n_spls_arr = vec_to_array (n_spls[i]);
        _delay.get (tailv, n_spls_arr);
        tap_tail[i] = vec1_array_unwrap (tailv);
        // feedback gain
        for (auto& val : tap_tail[i]) {
          val *= _param.fb_gain;
        }
      }
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
      case m_stereo2:
        stereo_2_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
      case m_ping_pong:
        pingpong_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
      case m_123:
      case m_132:
        x3_interleaving<T> (tap_head.data(), block, tap_tail.data());
        break;
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
      for (uint i = 0; i < block.size(); ++i) {
        auto ap_lfo_flt = _mod_lfo.tick_sine();
        ap_lfo_flt += 1; // unipolar downwards 32 samples
        ap_lfo_flt *= 0.5f * 32.f;
        auto ap_lfo = vec_to_array (vec_cast<uint> (ap_lfo_flt));

        for (uint j = 0; j < _diffusor.size(); ++j) {
          auto spl = make_vec (tap_head[i][j]);
          for (uint k = 0; k < _diffusor[0].size(); ++k) {
            spl = _diffusor[j][k].tick (
              spl,
              allpass_sizes[j][k] - ap_lfo[j],
              make_vec (_extpar.diffusion));
          }
          tap_head[i][j] = spl[0];
        }
      }
      // filter and feed back the samples
      for (uint i = 0; i < block.size(); ++i) {
        auto filt_x4   = vec_from_array (tap_head[i]);
        filt_x4        = _filters.tick<lp_idx> (filt_x4);
        auto head_vec1 = vec1_array_wrap (vec_to_array (filt_x4));
        _delay.push (head_vec1);
      }

      // pan
      std::array<std::array<float, 2>, n_taps> tap_mul;

      // specific output selection
      switch (_extpar.mode) {
      case m_stereo:
      case m_stereo2:
      case m_ping_pong:
        tap_mul[0] = get_pan (0.f);
        tap_mul[1] = get_pan (1.f);
        tap_mul[2] = get_pan (0.f);
        tap_mul[3] = get_pan (1.f);
        break;
      case m_123:
        tap_mul[0] = get_pan (0.f, 0.5f);
        tap_mul[1] = get_pan (0.f, 0.5f);
        tap_mul[2] = get_pan (0.5f);
        tap_mul[3] = get_pan (1.f);
        break;
      case m_132:
        tap_mul[0] = get_pan (0.f, 0.5f);
        tap_mul[1] = get_pan (0.f, 0.5f);
        tap_mul[2] = get_pan (1.f);
        tap_mul[3] = get_pan (0.5f);
        break;
      case m_1234:
        tap_mul[0] = get_pan (0.f);
        tap_mul[1] = get_pan (0.333333f);
        tap_mul[2] = get_pan (0.666666f);
        tap_mul[3] = get_pan (1.f);
        break;
      case m_2413:
        tap_mul[0] = get_pan (0.333333f);
        tap_mul[1] = get_pan (1.f);
        tap_mul[2] = get_pan (0.f);
        tap_mul[3] = get_pan (0.666666f);
        break;
      case m_1324:
        tap_mul[0] = get_pan (0.f);
        tap_mul[1] = get_pan (0.666666f);
        tap_mul[2] = get_pan (0.333333f);
        tap_mul[3] = get_pan (1.f);
        break;
      case m_2314:
        tap_mul[0] = get_pan (0.333333f);
        tap_mul[1] = get_pan (0.666666f);
        tap_mul[2] = get_pan (0.f);
        tap_mul[3] = get_pan (1.f);
        break;
      case m_1423:
        tap_mul[0] = get_pan (0.f);
        tap_mul[1] = get_pan (1.f);
        tap_mul[2] = get_pan (0.333333f);
        tap_mul[3] = get_pan (0.666666f);
        break;
      case m_chorus:
        tap_mul[0] = get_pan (0.f, 0.25f);
        tap_mul[1] = get_pan (0.333333f, 0.25f);
        tap_mul[2] = get_pan (0.666666f, 0.25f);
        tap_mul[3] = get_pan (1.f, 0.25f);
        break;
      default:
        assert (false);
        break;
      }
      // final tap acummulaton + gain
      for (uint i = 0; i < block.size(); ++i) {
        block[i][0] = 0.f;
        block[i][1] = 0.f;
        for (uint t = 0; t < n_taps; ++t) {
          block[i][0] += tap_tail[i][t] * tap_mul[t][0];
          block[i][1] += tap_tail[i][t] * tap_mul[t][1];
        }
        float gr = get_ducker_gain (block[i][0] + block[i][1]);
        block[i][0] *= _param.main_gain * gr;
        block[i][1] *= _param.main_gain * gr;
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
  static constexpr std::array<float, 2> get_pan (
    float pan,
    float correction = 1.f)
  {
    return {
      (float) (gcem::sin (M_PI_2 * (1.0f - pan)) * M_SQRT2 * correction),
      (float) (gcem::sin (M_PI_2 * pan) * M_SQRT2 * correction)};
  }
  //----------------------------------------------------------------------------
  template <class T>
  void stereo_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be inserted
    crange<std::array<T, 2> const>  in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    // stereo with side ping pong
    auto angle = get_fdn4_angle();
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][0];
      tap_head[i][1] = tap_tail[i][1];
      tap_head[i][2] = tap_tail[i][2];
      tap_head[i][3] = tap_tail[i][3];
      // Extra FDN diffusion. Mostly will be noticed with allpass modulation
      tap_head[i] = rotation_matrix<4>::tick<arith_type> (tap_head[i], angle);
      tap_head[i][0] += in[i][0];
      tap_head[i][2] += in[i][1];
      auto side = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][1] += side;
      tap_head[i][3] += side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void stereo_2_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be inserted
    crange<std::array<T, 2> const>  in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    // stereo with side ping pong
    auto angle = get_fdn4_angle();
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][0];
      tap_head[i][1] = tap_tail[i][1];
      tap_head[i][2] = tap_tail[i][3];
      tap_head[i][3] = tap_tail[i][2];
      // Extra FDN diffusion. Mostly will be noticed with allpass modulation
      tap_head[i] = rotation_matrix<4>::tick<arith_type> (tap_head[i], angle);
      tap_head[i][0] += in[i][0];
      tap_head[i][2] += in[i][1];
      auto side = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][1] += side;
      tap_head[i][3] += side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void pingpong_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be inserted
    crange<std::array<T, 2> const>  in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][1];
      tap_head[i][1] = tap_tail[i][0];
      tap_head[i][2] = tap_tail[i][3];
      tap_head[i][3] = tap_tail[i][2];
      auto mid       = (arith_type) (in[i][0] + in[i][1]);
      auto side      = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][0] += mid;
      tap_head[i][2] += side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void x3_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be inserted
    crange<std::array<T, 2> const>  in, // inputs
    std::array<arith_type, n_taps> const* tap_tail) // samples that are outputs
  {
    for (uint i = 0; i < in.size(); ++i) {
      tap_head[i][0] = tap_tail[i][3] * 0.5f;
      tap_head[i][1] = tap_tail[i][3] * 0.5f;
      tap_head[i][2] = tap_tail[i][0] + tap_tail[i][1];
      tap_head[i][3] = tap_tail[i][2];
      auto mid       = (arith_type) (in[i][0] + in[i][1]);
      auto side      = (arith_type) (in[i][0] - in[i][1]);
      tap_head[i][0] += mid;
      tap_head[i][1] += side;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void x4_interleaving (
    std::array<arith_type, n_taps>* tap_head, // samples that will be inserted
    crange<std::array<T, 2> const>  in, // inputs
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
    std::array<arith_type, n_taps>* tap_head, // samples that will be inserted
    crange<std::array<T, 2> const>  in, // inputs
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
  struct external_parameters {
    uint  mode;
    uint  mod_mode;
    float delay_spls;
    float feedback;
    float diffusion;
    float gain;
    float ducking_threshold;
    float ducking_speed;
    float mode_param;
    float mod_freq;
    float mod_depth;
    float mod_spread;
    float damp_freq;
    float tilt_db;
    float stereo;
  };
  //----------------------------------------------------------------------------
  struct internal_parameters {
    double delay_spls_max;
    double spls_x_beat;
    float  fb_gain;
    float  main_gain;
  };
  //----------------------------------------------------------------------------
  std::array<std::array<allpass<float_x1>, 4>, n_taps> _diffusor;
  //----------------------------------------------------------------------------
  external_parameters                           _extpar {};
  internal_parameters                           _param {};
  block_resampler<arith_type, 2>                _resampler {};
  modulable_thiran2_delay_line<vec1_type>       _delay {};
  std::vector<vec1_type>                        _mem {};
  part_class_array<onepole_smoother, vec1_type> _n_spls_smoother {};
  part_class_array<tilt_eq, double_x2>          _tilt {};
  lfo<n_taps>                                   _mod_lfo;
  lfo<n_taps>                                   _ap_lfo;
  part_class_array<slew_limiter, vec1_type>     _ducker_follow;
  enum { dc_idx, lp_idx };
  part_classes<
    mp_list<mystran_dc_blocker, onepole_lowpass>,
    vec<arith_type, n_taps>,
    false>
    _filters;
};
//------------------------------------------------------------------------------
} // namespace artv
