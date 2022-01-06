#pragma once

#include <complex>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/biquad.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros_reversed.hpp"

#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
class lin_eq4x {
private:
  enum class paramtype { band_type, frequency, q, gain, diff, quality };

public:
  template <uint Band, paramtype Type>
  struct param {
    static constexpr uint band  = Band;
    static constexpr uint ptype = (uint) Type;
  };
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::eq;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    for (auto& eq : _eq) {
      eq.reset(); // realloc
    }
    using x1_t = vec<decltype (_smooth_coeff), 1>;
    _sample    = 0;
    _latency   = 0;
    _quality   = 0;

    _plugcontext = &pc;
    _cfg         = decltype (_cfg) {};

    reset_quality_settings (pc.get_sample_rate());

    smoother::reset_coeffs (
      make_crange (_smooth_coeff).cast (x1_t {}),
      vec_set<x1_t> (1 / 0.001),
      pc.get_sample_rate());
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    uint n_enabled_bands = 0;

    std::array<bool, 2>       is_copied {false, false};
    std::array<double_x2, 16> io;

    for (uint b = 0; b < n_bands; ++b) {
      if (_cfg[b].has_changes) {
        reset_band (b);
      }
      if (_cfg[b].type == bandtype::off) {
        continue;
      }

      ++n_enabled_bands;

      for (uint offset = 0; offset < samples; offset += io.size()) {
        uint blocksize = std::min<uint> (io.size(), samples - offset);
        // interleaving

        std::array<T const*, 2> src;
        src[0] = is_copied[0] ? &outs[0][offset] : &ins[0][offset];
        src[1] = is_copied[1] ? &outs[1][offset] : &ins[1][offset];

        for (uint i = 0; i < blocksize; ++i) {
          io[i][0] = src[0][i];
          io[i][1] = src[1][i];
        }
        // bulk smoothing (at blocksize rate)
        crange<double_x2> internal = _eq[b].get_all_coeffs();
        for (uint j = 0; j < _target_coeffs[b].size(); ++j) {
          for (uint i = 0; i < blocksize; ++i) {
            internal[j] = smoother::tick (
              make_crange (_smooth_coeff),
              make_crange (internal[j]),
              _target_coeffs[b][j]);
          }
        }

        // backward zeros
        for (uint i = 0; i < blocksize; ++i) {
          io[i]
            = _eq[b].tick<bckwd_zeros_idx> (io[i], biquad::rev_zeros_tag {});
        }
        // backward poles
        _eq[b].tick<bckwd_poles_idx> (
          make_crange (io.data(), blocksize),
          _n_stages[b][_quality],
          _sample + offset);
        // forward
        for (uint i = 0; i < blocksize; ++i) {
          io[i] = _eq[b].tick<fwd_idx> (io[i]);
        }
        // deinterleaving
        switch (_topology) {
        case topology::stereo:
          for (uint i = 0; i < blocksize; ++i) {
            outs[0][offset + i] = io[i][0];
            outs[1][offset + i] = io[i][1];
          }
          break;
        case topology::l:
          for (uint i = 0; i < blocksize; ++i) {
            outs[0][offset + i] = io[i][0];
          }
          break;
        case topology::r:
          for (uint i = 0; i < blocksize; ++i) {
            outs[1][offset + i] = io[i][1];
          }
          break;
        default:
          assert (false);
          break;
        }
      }

      switch (_topology) {
      case topology::stereo:
        is_copied[0] = is_copied[1] = true;
        break;
      case topology::l:
        is_copied[0] = true;
        break;
      case topology::r:
        is_copied[1] = true;
        break;
      default:
        break;
      }
    }
    for (uint c = 0; c < 2; ++c) {
      if (unlikely (!is_copied[c])) {
        memcpy (outs[c], ins[c], samples * sizeof outs[c][0]);
      }
    }
    _sample += samples;
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::band_type>, int v)
  {
    static_assert (band < n_bands, "");
    if (v >= (int) bandtype::size || v < 0) {
      v = (int) bandtype::off;
    }
    if ((int) _cfg[band].type == (int) v) {
      return;
    }
    _cfg[band].has_changes      = true;
    _cfg[band].reset_band_state = true;
    _cfg[band].type             = (bandtype) v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::band_type>)
  {
    return choice_param (
      0, make_cstr_array ("Off", "Peak", "LowShelf", "HighShelf"), 16);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::frequency>, float v)
  {
    static_assert (band < n_bands, "");
    if (_cfg[band].freq_note == v) {
      return;
    }
    _cfg[band].has_changes = true;
    // the filter generates spikes on its delay line if not sweeped gently,
    // protect from them by doing a full clear on abrupt note sweeps. Notice
    // that it's only sweeping from high to low frequency that is a problem.
    _cfg[band].has_changes |= (_cfg[band].freq_note - v) > 12.;
    _cfg[band].freq_note = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::frequency>)
  {
    if constexpr (band == 0) {
      return frequency_parameter (80.0, 20000.0, 100.0);
    }
    else if constexpr (band == 1) {
      return frequency_parameter (200.0, 20000.0, 300.0);
    }
    else if constexpr (band == 2) {
      return frequency_parameter (600.0, 20000.0, 1000.0);
    }
    else if constexpr (band == 3) {
      return frequency_parameter (2000.0, 20000.0, 3000.0);
    }
    else {
      return frequency_parameter (20.0, 20000.0, 440.0);
    }
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::q>, float v)
  {
    static_assert (band < n_bands, "");

    _cfg[band].has_changes |= _cfg[band].q != v;
    _cfg[band].q = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::q>)
  {
    if constexpr (band == 0) {
      return float_param ("", 0.1f, 1.0f, 0.5f, 0.0001f, 0.9f);
    }
    else if constexpr (band == 1) {
      return float_param ("", 0.1f, 1.7f, 0.5f, 0.0001f, 0.8f);
    }
    else if constexpr (band == 2) {
      return float_param ("", 0.1f, 3.f, 0.5f, 0.0001f, 0.6f);
    }
    else if constexpr (band == 3) {
      return float_param ("", 0.1f, 6.f, 0.5f, 0.0001f, 0.3f);
    }
    else {
      return float_param ("", 0.1f, 2.f, 0.5f, 0.0001f, 0.3f);
    }
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::gain>, float v)
  {
    static_assert (band < n_bands, "");
    v *= 0.5; // cascaded response...
    _cfg[band].has_changes |= _cfg[band].gain_db != v;
    _cfg[band].gain_db = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::gain>)
  {
    return float_param ("dB", -20.f, 20.f, 0.f, 0.25f, 0.6f, true);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::diff>, float v)
  {
    static_assert (band < n_bands, "");
    v *= 0.01 * 0.25; // a quarter octave up and down
    _cfg[band].has_changes |= _cfg[band].diff != v;
    _cfg[band].diff = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::diff>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.05f);
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_quality_steps = 4;

  struct quality_tag {};

  void set (quality_tag, int v)
  {
    // Can't happen while processing!
    assert (v >= 0 && v < n_quality_steps);
    if (_quality == v) {
      return;
    }
    for (auto& band : _cfg) {
      band.has_changes      = true;
      band.reset_band_state = true;
    }
    _quality = v;
  }

  static constexpr auto get_parameter (quality_tag)
  {
    return choice_param (
      0, make_cstr_array ("Normal", "High", "Overspeced", "Snake oil"), 8);
  }
  //----------------------------------------------------------------------------
  struct topology_tag {};
  void set (topology_tag, int v) { _topology = (topology) v; }

  static constexpr auto get_parameter (topology_tag)
  {
    return choice_param (0, make_cstr_array ("Stereo", "L", "R"), 16);
  }
  //----------------------------------------------------------------------------
  using band1_type_tag = param<0, paramtype::band_type>;
  using band2_type_tag = param<1, paramtype::band_type>;
  using band3_type_tag = param<2, paramtype::band_type>;
  using band4_type_tag = param<3, paramtype::band_type>;
  using band1_freq_tag = param<0, paramtype::frequency>;
  using band2_freq_tag = param<1, paramtype::frequency>;
  using band3_freq_tag = param<2, paramtype::frequency>;
  using band4_freq_tag = param<3, paramtype::frequency>;
  using band1_q_tag    = param<0, paramtype::q>;
  using band2_q_tag    = param<1, paramtype::q>;
  using band3_q_tag    = param<2, paramtype::q>;
  using band4_q_tag    = param<3, paramtype::q>;
  using band1_gain_tag = param<0, paramtype::gain>;
  using band2_gain_tag = param<1, paramtype::gain>;
  using band3_gain_tag = param<2, paramtype::gain>;
  using band4_gain_tag = param<3, paramtype::gain>;
  using band1_diff_tag = param<0, paramtype::diff>;
  using band2_diff_tag = param<1, paramtype::diff>;
  using band3_diff_tag = param<2, paramtype::diff>;
  using band4_diff_tag = param<3, paramtype::diff>;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    band1_type_tag,
    band1_freq_tag,
    band1_gain_tag,
    band1_q_tag,
    band2_type_tag,
    band2_freq_tag,
    band2_gain_tag,
    band2_q_tag,
    band3_type_tag,
    band3_freq_tag,
    band3_gain_tag,
    band3_q_tag,
    band4_type_tag,
    band4_freq_tag,
    band4_gain_tag,
    band4_q_tag,
    band1_diff_tag,
    band2_diff_tag,
    band3_diff_tag,
    band4_diff_tag,
    quality_tag,
    topology_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void try_reset_latency()
  {
    uint new_latency = 0;
    for (uint b = 0; b < n_bands; ++b) {
      uint latency = (1 << _n_stages[b][_quality]) + 1;
      new_latency += _cfg[b].type == bandtype::off ? 0 : latency;
    }
    if (_latency != new_latency) {
      _latency = new_latency;
      _plugcontext->set_delay_compensation (new_latency);
    }
  }
  //----------------------------------------------------------------------------
  enum class bandtype { off, svf_bell, svf_lshelf, svf_hshelf, size };
  //----------------------------------------------------------------------------
  struct bandconfig {
    bandtype type             = bandtype::off;
    float    freq_note        = constexpr_midi_note_to_hz (440.f);
    uint     freq_spl         = 0;
    float    q                = 0.7f;
    float    gain_db          = 0.f;
    float    diff             = 0.f;
    bool     has_changes      = false;
    bool     reset_band_state = false;
  };
  //----------------------------------------------------------------------------
  void reset_band (uint band)
  {
    auto& b             = _cfg[band];
    auto  bandtype_prev = b.type;
    auto  sr            = (float) _plugcontext->get_sample_rate();

    auto freq = vec_set<double_x2> (midi_note_to_hz (b.freq_note));
    freq[1] *= exp2 (b.diff);
    // not changing the Q, as we want both poles to either be real or conjugate
    auto q    = vec_set<double_x2> (b.q);
    auto gain = vec_set<double_x2> (b.gain_db);

    auto co        = make_crange (_target_coeffs[band]);
    auto co_biquad = co.shrink_head (_eq[0].get_coeff_offset<fwd_idx>());
    auto co_bwd_poles
      = co.shrink_head (_eq[0].get_coeff_offset<bckwd_poles_idx>());
    auto co_bwd_zeros
      = co.shrink_head (_eq[0].get_coeff_offset<bckwd_zeros_idx>());

    std::array<double_x2, biquad::n_coeffs> coeffs;

    switch (b.type) {
    case bandtype::off:
      _eq[band].reset_states<fwd_idx>();
      _eq[band].reset_states<bckwd_poles_idx> (_n_stages[band][_quality]);
      _eq[band].reset_states<bckwd_zeros_idx>();
      break;
    case bandtype::svf_bell:
      _eq[band].reset_coeffs_ext<fwd_idx> (
        co_biquad, freq, q, gain, sr, bell_tag {});
      break;
    case bandtype::svf_lshelf:
      _eq[band].reset_coeffs_ext<fwd_idx> (
        co_biquad, freq, q, gain, sr, lowshelf_tag {});
      break;
    case bandtype::svf_hshelf:
      _eq[band].reset_coeffs_ext<fwd_idx> (
        co_biquad, freq, q, gain, sr, highshelf_tag {});
      break;
    default:
      jassert (false);
    }

    crange_memcpy (co_bwd_zeros, co_biquad);

    std::array<vec_complex<double_x2>, 2> poles;
    biquad::get_poles<double_x2> (co_biquad, poles);

    _eq[band].reset_coeffs_ext<bckwd_poles_idx> (
      co_bwd_poles, poles[0], poles[1]);

    // reset smoothing
    if (_cfg[band].reset_band_state) {
      _eq[band].reset_states<bckwd_poles_idx> (_n_stages[band][_quality]);
      _eq[band].reset_states<bckwd_zeros_idx>();
      _eq[band].reset_states<fwd_idx>();
      crange_copy<double_x2> (_eq[band].get_all_coeffs(), _target_coeffs[band]);
      try_reset_latency();
    }
    _cfg[band].has_changes      = false;
    _cfg[band].reset_band_state = false;
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_bands    = 4;
  static constexpr uint n_channels = 2;
  //----------------------------------------------------------------------------
  static uint get_n_stages (
    double freq,
    double q,
    double gain,
    double samplerate,
    double snr_db)
  {
    std::array<double_x1, biquad::n_coeffs> co;
    biquad::reset_coeffs<double_x1> (
      co,
      vec_set<1> (freq),
      vec_set<1> (q),
      vec_set<1> (gain),
      samplerate,
      bell_tag {});
    std::array<vec_complex<double_x1>, 2> poles;
    biquad::get_poles<double_x1> (co, poles);
    return get_reversed_pole_n_stages (poles[0].to_std (0), snr_db);
  }
  //----------------------------------------------------------------------------
  void reset_quality_settings (double samplerate)
  {
    std::array<double, n_bands> freq, gain, q;

    mp11::mp_for_each<mp11::mp_iota_c<n_bands>> ([&] (auto band) {
      constexpr uint b = decltype (band)::value;

      using freq_param = param<b, paramtype::frequency>;
      using q_param    = param<b, paramtype::q>;
      using gain_param = param<b, paramtype::gain>;

      freq[b] = get_parameter (freq_param {}).min;
      freq[b] = midi_note_to_hz (freq[b]);
      gain[b] = get_parameter (gain_param {}).max;
      q[b]    = get_parameter (q_param {}).max;
    });

    for (uint b = 0; b < n_bands; ++b) {
      _n_stages[b][0] = get_n_stages (freq[b], q[b], gain[b], samplerate, 90.);
      _n_stages[b][1] = get_n_stages (freq[b], q[b], gain[b], samplerate, 120.);
      _n_stages[b][2] = _n_stages[b][1] + 1;
      _n_stages[b][3] = _n_stages[b][1] + 2;
      static_assert (n_quality_steps == 4, "Update this!");

      for (auto& stages : _n_stages[b]) {
        stages = std::min (stages, max_n_stages);
      }
    }
  }
  //----------------------------------------------------------------------------
  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  enum class topology { stereo, l, r };
  //----------------------------------------------------------------------------
  // Big amount of preallocated memory to cope with high samplerates. This was
  // to use the convenience of "part_classes" at the expense of always
  // allocating the worst case amount of memory.
  static constexpr uint max_n_stages = 18;
  //----------------------------------------------------------------------------
  enum idx {
    bckwd_poles_idx,
    bckwd_zeros_idx,
    fwd_idx,
  };

  using eqs = part_classes<
    mp_list<
      // make_max_stages_t_rev<t_rev_cpole_pair_czero_pair, max_n_stages>,
      make_max_stages_t_rev<t_rev_pole_pair, max_n_stages>,
      biquad,
      biquad>,

    double_x2,
    true>;

  using coeff_array = std::array<eqs::value_type, eqs::n_coeffs>;

  alignas (eqs::value_type) std::array<coeff_array, n_bands> _target_coeffs;
  std::array<eqs, n_bands> _eq;

  std::array<bandconfig, n_bands>                        _cfg;
  uint                                                   _sample;
  double                                                 _smooth_coeff;
  topology                                               _topology;
  uint                                                   _quality;
  std::array<std::array<uint, n_quality_steps>, n_bands> _n_stages;
  uint                                                   _latency     = 0;
  plugin_context*                                        _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
