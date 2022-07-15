#pragma once

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
    v *= v;
    v *= 0.5f; // cut at 0.5f for the allpass
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
  void set (ducking_threshold_tag, float v)
  {
    // TODO: dB!
    if (v == _extpar.ducking_threshold) {
      return;
    }
    _extpar.ducking_threshold = v;
  }

  static constexpr auto get_parameter (ducking_threshold_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
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
    // TODO:
  }

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Ping-Pong",
        "Pong-Ping",
        "Stereo",
        "L-C-R",
        "R-C-L",
        "L-2C-R",
        "R-2C-L"
        "Rotating"),
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
    _mod_lfo.set_freq (vec_set<n_delay_lines> (v), (float) sr_target_freq);
  }

  static constexpr auto get_parameter (mod_freq_tag)
  {
    return float_param ("Hz", 0.f, 10.f, 0.f, 0.01f, 0.5f);
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
  struct mod_spread_tag {};
  void set (mod_spread_tag, float v)
  {
    v *= 0.01f;
    if (v == _extpar.mod_spread) {
      return;
    }
    _extpar.mod_spread = v;
  }

  static constexpr auto get_parameter (mod_spread_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
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
      0, make_cstr_array ("Random", "Chorus1", "Chorus2", "Chorus3"), 10);
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
    return float_param ("dB", -12.f, 12.f, 0.f, 0.1f);
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
    ducking_speed_tag,
    ducking_threshold_tag,
    mod_freq_tag,
    mod_depth_tag,
    mod_mode_tag,
    mod_spread_tag>;
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
      vec_set<n_delay_lines> (1.f), (float) sr_target_freq);

    _filters.reset_states<lp_idx>();
    _filters.reset_states<dc_idx>();

    _tilt.reset_states();

    using phase_type = decltype (_mod_lfo)::phase_type;
    using value_type = decltype (_mod_lfo)::value_type;

    _mod_lfo.reset();
    _mod_lfo.set_phase (phase_type {
      value_type {0.f, 0.25f, 0.5f, 0.75f}, phase_type::normalized {}});

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
  static constexpr uint blocksize       = 32;
  static constexpr uint n_delay_lines   = 4;
  static constexpr uint max_mod_samples = 132;
  using arith_type                      = float;
  using vec1_type                       = vec<float, 1>;
  //----------------------------------------------------------------------------
  void initialize_buffer_related_parts()
  {
    // compute memory requirements
    double max_spls = max_mod_samples;
    max_spls += _param.spls_x_beat * max_t_beats;
    auto fdn_size         = pow2_round_ceil ((uint) std::ceil (max_spls));
    _param.delay_spls_max = (double) (fdn_size);
    uint n_samples_delay  = _delay.n_required_elems (n_delay_lines, fdn_size);

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
    _delay.reset (mem.cut_head (n_samples_delay), n_delay_lines);

    for (uint i = 0; i < allpass_sizes.size(); ++i) {
      auto& line_sizes = allpass_sizes[i];
      for (uint j = 0; j < line_sizes.size(); ++j) {
        _diffusor[i][j].reset (mem.cut_head (allpass_sizes[i][j]));
      }
    }
  }
  //----------------------------------------------------------------------------
  static std::array<std::array<uint, 4>, n_delay_lines>
  get_diffusor_delay_spls()
  {
    return {
      {{{225u, 556u, 441u, 341u}},
       {{225u, 556u, 441u, 341u}},
       {{225u, 556u, 441u, 341u}},
       {{225u, 556u, 441u, 341u}}}};
  }
  //----------------------------------------------------------------------------
  static std::array<uint, n_delay_lines> get_diffusor_delay_total_spls()
  {
    auto                            spls = get_diffusor_delay_spls();
    std::array<uint, n_delay_lines> ret;
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
  template <class T>
  void process_block (crange<std::array<T, 2>> io)
  {
    // TODO: when do, split in sublocks for performance's sake.
    while (io.size()) {
      auto block = io.cut_head (std::min<uint> (io.size(), blocksize));
      for (uint i = 0; i < block.size(); ++i) {

        auto n_spls_smooth
          = _n_spls_smoother.tick (vec_set<1> (_extpar.delay_spls))[0];
        vec<arith_type, n_delay_lines> n_spls
          = vec_set<n_delay_lines> ((arith_type) n_spls_smooth);

        auto diffusor_correction = vec_cast<arith_type> (
          vec_from_array (get_diffusor_delay_total_spls()));

        n_spls -= diffusor_correction;

        // delay lfo
        vec<arith_type, n_delay_lines> mod;
        switch (_extpar.mod_mode) {
        case 0:
          mod = _mod_lfo.tick_filt_sample_and_hold();
          break;
        case 1:
          mod = _mod_lfo.tick_sine();
          break;
        case 2:
          mod = _mod_lfo.tick_triangle();
          break;
        case 3:
          mod = _mod_lfo.tick_trapezoid (vec_set<n_delay_lines> (0.75f));
          break;
        default:
          assert (false);
          break;
        };
        mod *= (arith_type) (max_mod_samples * _extpar.mod_depth);
        n_spls += mod;
        n_spls = vec_min (2.f, n_spls); // 2 samples for the thiran2 interp

        auto& sample = block[i];

        // tilt
        double_x2 tilted {sample[0], sample[1]};
        tilted    = _tilt.tick (tilted);
        sample[0] = tilted[0];
        sample[1] = tilted[1];

        // ping pong.
        auto mid  = (arith_type) (sample[0] + sample[1]);
        auto side = (arith_type) (sample[0] - sample[1]);

        std::array<arith_type, n_delay_lines> spls_tail;
        std::array<vec1_type, n_delay_lines>  spls_tail_vec1;
        auto n_spls_arr = vec_to_array (n_spls);
        _delay.get (spls_tail_vec1, n_spls_arr);
        spls_tail = vec1_array_unwrap (spls_tail_vec1);

        std::array<arith_type, n_delay_lines> spls_head;
        spls_head[0] = spls_tail[1];
        spls_head[1] = spls_tail[0];
        spls_head[2] = spls_tail[3];
        spls_head[3] = spls_tail[2];

        for (auto& v : spls_head) {
          v *= _param.fb_gain;
        }

        spls_head[0] += mid;
        spls_head[2] += side;

        // difussion
        auto allpass_sizes = get_diffusor_delay_spls();
        for (uint j = 0; j < _diffusor.size(); ++j) {
          auto diffuse = make_vec (spls_head[j]);
          for (uint k = 0; k < _diffusor[0].size(); ++k) {
            diffuse = _diffusor[j][k].tick (
              diffuse, allpass_sizes[j][k], make_vec (_extpar.diffusion));
          }
          spls_head[j] = diffuse[0];
        }

        // filters
        auto filt_x4 = vec_from_array (spls_head);
        filt_x4      = _filters.tick<lp_idx> (filt_x4);
        filt_x4      = _filters.tick<dc_idx> (filt_x4);

        // insert
        spls_head           = vec_to_array (filt_x4);
        auto spls_head_vec1 = vec1_array_wrap (spls_head);
        _delay.push (spls_head_vec1);

        sample[0] = spls_tail[0] + spls_tail[2];
        sample[1] = spls_tail[1] + spls_tail[3];

        sample[0] *= _param.main_gain;
        sample[1] *= _param.main_gain;
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
  // starting point freeverb diffusor
  std::array<std::array<allpass<float_x1>, 4>, n_delay_lines> _diffusor;
  //----------------------------------------------------------------------------
  external_parameters                           _extpar {};
  internal_parameters                           _param {};
  block_resampler<arith_type, 2>                _resampler {};
  modulable_thiran2_delay_line<vec1_type>       _delay {};
  std::vector<vec1_type>                        _mem {};
  part_class_array<onepole_smoother, vec1_type> _n_spls_smoother {};
  part_class_array<tilt_eq, double_x2>          _tilt {};
  lfo<n_delay_lines>                            _mod_lfo;
  enum { dc_idx, lp_idx };
  part_classes<
    mp_list<mystran_dc_blocker, onepole_lowpass>,
    vec<arith_type, n_delay_lines>,
    false>
    _filters;
};
//------------------------------------------------------------------------------
} // namespace artv
