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
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

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
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    bool in_on_out = false;

    std::array<bool, 2>       is_copied {false, false};
    std::array<double_x2, 32> io;
    std::array<double_x2, 32> in_cp;
    std::array<double_x2, 32> sum;

    // recompute coefficients
    for (uint b = 0; b < n_bands; ++b) {
      if (_cfg[b].has_changes) {
        reset_band (b);
      }
    }
    // process serial bands
    for (uint bi = 0; (bi < n_bands) && (_order_serial[bi] >= 0); ++bi) {
      auto b = _order_serial[bi];
      for (uint offset = 0; offset < samples; offset += io.size()) {
        uint blocksize = std::min<uint> (io.size(), samples - offset);

        // interleaving
        std::array<T const*, 2> src;
        src[0] = in_on_out ? &outs[0][offset] : &ins[0][offset];
        src[1] = in_on_out ? &outs[1][offset] : &ins[1][offset];

        for (uint i = 0; i < blocksize; ++i) {
          io[i][0] = src[0][i];
          io[i][1] = src[1][i];
        }
        process_band (make_crange (io.data(), blocksize), b, _sample + offset);
        // deinterleaving
        for (uint i = 0; i < blocksize; ++i) {
          outs[0][offset + i] = io[i][0];
          outs[1][offset + i] = io[i][1];
        }
      }
      in_on_out = true;
    }
    // early exit if no parallel bands
    if (_order_parallel[0] < 0) {
      if (unlikely (!in_on_out)) {
        memcpy (outs[0], ins[0], samples * sizeof outs[0][0]);
        memcpy (outs[1], ins[1], samples * sizeof outs[1][0]);
      }
      _sample += samples;
      return;
    }
    // process parallel bands
    for (uint offset = 0; offset < samples; offset += io.size()) {
      uint blocksize = std::min<uint> (io.size(), samples - offset);

      std::array<T const*, 2> src;
      src[0] = in_on_out ? &outs[0][offset] : &ins[0][offset];
      src[1] = in_on_out ? &outs[1][offset] : &ins[1][offset];

      // interleaving
      for (uint i = 0; i < blocksize; ++i) {
        in_cp[i][0] = src[0][i];
        in_cp[i][1] = src[1][i];
      }

      // do a parallel sum of all the bands
      for (uint bi = 0; (bi < n_bands) && (_order_parallel[bi] >= 0); ++bi) {
        auto b = _order_parallel[bi];

        auto buff = make_crange ((bi == 0) ? &sum[0] : &io[0], blocksize);
        crange_memcpy (buff, make_crange (&in_cp[0], blocksize));
        process_band (buff, b, _sample + offset);

        // match latency with other bands
        if (_bandcomp[b].delay() != 0) {
          for (uint i = 0; i < blocksize; ++i) {
            buff[i] = _bandcomp[b].exchange (buff[i]);
          }
        }
        // scale
        for (uint i = 0; i < blocksize; ++i) {
          buff[i] *= _cfg[b].slope_gain;
        }
        // append results
        if (bi != 0) {
          for (uint i = 0; i < blocksize; ++i) {
            sum[i] += buff[i];
          }
        }
      }
      // sum with input
      for (uint i = 0; i < blocksize; ++i) {
        in_cp[i] = _incomp.exchange (in_cp[i]);
        in_cp[i] += sum[i];
      }
      // deinterleave
      for (uint i = 0; i < blocksize; ++i) {
        outs[0][offset + i] = in_cp[i][0];
        outs[1][offset + i] = in_cp[i][1];
      }
    }
    _sample += samples;
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::band_type>, int v)
  {
    static_assert (band < n_bands, "");
    if (v >= (int) bandtype::count || v < 0) {
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
      0,
      make_cstr_array (
        "Off", "Peak", "LowShelf", "HighShelf", "LowPass", "HighPass"),
      16);
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
      return frequency_parameter (470.0, 20000.0, 1000.0);
    }
    else if constexpr (band == 3) {
      return frequency_parameter (800.0, 20000.0, 3000.0);
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
      return float_param ("", 0.1f, 1.0f, M_SQRT1_2, 0.0001f, 0.9f);
    }
    else if constexpr (band == 1) {
      return float_param ("", 0.1f, 2.5f, M_SQRT1_2, 0.0001f, 0.8f);
    }
    else if constexpr (band == 2) {
      return float_param ("", 0.1f, 3.f, M_SQRT1_2, 0.0001f, 0.6f);
    }
    else if constexpr (band == 3) {
      return float_param ("", 0.1f, 5.f, M_SQRT1_2, 0.0001f, 0.3f);
    }
    else {
      return float_param ("", 0.1f, 2.f, M_SQRT1_2, 0.0001f, 0.3f);
    }
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::gain>, float v)
  {
    static_assert (band < n_bands, "");
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
  static constexpr uint n_quality_steps = 3;

  struct quality_tag {};

  void set (quality_tag, int v)
  {
    // Can't happen while processing!
    if (_quality == v || v >= n_quality_steps) {
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
      0, make_cstr_array ("Normal", "Very-High", "Snake oil"), 8);
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
    quality_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    for (auto& eq : _eq) {
      eq.reset(); // realloc
    }
    clear_orders();

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

    mp11::mp_for_each<parameters> ([=] (auto param) {
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void process_band (crange<double_x2> io, uint b, uint sample_idx)
  {
    // bulk smoothing (at blocksize rate)
    crange<double_x2> internal = _eq[b].get_all_coeffs();
    for (uint j = 0; j < _target_coeffs[b].size(); ++j) {
      for (uint i = 0; i < io.size(); ++i) {
        internal[j] = smoother::tick (
          make_crange (_smooth_coeff),
          make_crange (internal[j]),
          _target_coeffs[b][j]);
      }
    }
    // backward zeros
    for (uint i = 0; i < io.size(); ++i) {
      io[i] = _eq[b].tick<bckwd_zeros_idx> (io[i], biquad::rev_zeros_tag {});
    }
    // backward poles
    _eq[b].tick<bckwd_poles_idx> (io, _n_stages[b][_quality], sample_idx);
    // forward
    for (uint i = 0; i < io.size(); ++i) {
      io[i] = _eq[b].tick<fwd_idx> (io[i]);
    }
  }
  //----------------------------------------------------------------------------
  void reset_latency()
  {
    // To minimize latency the EQ processes the shelves and peaks in parallel,
    // as it is mostly equivalent. This is not so with low and highpasses, so
    // the low and highpasses are processed first and they add to the latency
    // amount.
    clear_orders();
    uint order_idx   = 0;
    uint new_latency = 0;
    // get per-band latencies and all the serial bands
    for (uint b = 0; b < n_bands; ++b) {
      if (_cfg[b].type == bandtype::off) {
        _cfg[b].latency = 0;
        continue;
      }
      uint latency    = (1 << _n_stages[b][_quality]) + 1;
      _cfg[b].latency = latency;

      if (!is_parallel (_cfg[b].type)) {
        _order_serial[order_idx] = b;
        ++order_idx;
        new_latency += latency;
      }
    }

    order_idx             = 0;
    uint n_parallel       = 0;
    uint latency_parallel = 0;
    uint mem_size         = 0;

    // get the parallel bands
    for (uint b = 0; b < n_bands; ++b) {
      if (_cfg[b].type == bandtype::off) {
        continue;
      }
      if (is_parallel (_cfg[b].type)) {
        ++n_parallel;
        latency_parallel = std::max (_cfg[b].latency, latency_parallel);
        mem_size += (latency_parallel - _cfg[b].latency);
        _order_parallel[order_idx] = b;
        ++order_idx;
      }
    }
    mem_size += latency_parallel; // for the in buffer of the parallel section
    new_latency += latency_parallel;
    // adjust compensation and rebuild latency buffers
    _latency = new_latency;
    _plugcontext->set_delay_compensation (new_latency);
    _mem.clear();
    _mem.resize (mem_size);

    // assign the buffers to the parallel bands
    auto ptr = _mem.data();
    _incomp.reset (make_crange (ptr, latency_parallel), latency_parallel);
    ptr += latency_parallel;

    for (uint b = 0; b < n_bands; ++b) {
      if (!is_parallel (_cfg[b].type) || _cfg[b].type == bandtype::off) {
        continue; // series. no compensation required
      }
      uint comp = latency_parallel - _cfg[b].latency;
      _bandcomp[b].reset (make_crange (ptr, comp), comp);
      ptr += comp;
    }
  }
  //----------------------------------------------------------------------------
  void reset_band (uint band)
  {
    auto& b             = _cfg[band];
    auto  bandtype_prev = b.type;
    auto  sr            = (float) _plugcontext->get_sample_rate();

    auto freq = vec_set<double_x2> (midi_note_to_hz (b.freq_note));
    freq[1] *= exp2 (b.diff);
    // not changing the Q, as we want both poles to either be real or
    // conjugate
    auto q = vec_set<double_x2> (b.q);

    auto gain = exp ((double) b.gain_db * (1. / 20.) * M_LN10);

    b.slope_gain = gain - 1; // e.g 0dB = 1, so the multiplier is 0

    auto co        = make_crange (_target_coeffs[band]);
    auto co_biquad = co.advanced (_eq[0].get_coeff_offset<fwd_idx>());
    auto co_bwd_poles
      = co.advanced (_eq[0].get_coeff_offset<bckwd_poles_idx>());
    auto co_bwd_zeros
      = co.advanced (_eq[0].get_coeff_offset<bckwd_zeros_idx>());

    std::array<double_x2, biquad::n_coeffs> coeffs;
    bool high_srate = _plugcontext->get_sample_rate() > 70000;

    switch (b.type) {
    case bandtype::off:
      break;
    case bandtype::bell:
      q *= sqrt (gain);
      if (high_srate) {
        _eq[band].reset_coeffs_ext<fwd_idx> (
          co_biquad, freq, q, sr, bandpass_tag {});
      }
      else {
        _eq[band].reset_coeffs_ext<fwd_idx> (
          co_biquad, freq, q, sr, biquad::mvic_bandpass_hq_tag {});
      }
      break;
    case bandtype::lshelf:
    case bandtype::lp:
      if (high_srate) {
        _eq[band].reset_coeffs_ext<fwd_idx> (
          co_biquad, freq, q, sr, lowpass_tag {});
      }
      else {
        _eq[band].reset_coeffs_ext<fwd_idx> (
          co_biquad, freq, q, sr, biquad::mvic_lowpass_hq_tag {});
      }
      break;
    case bandtype::hshelf:
    case bandtype::hp:
      if (high_srate) {
        _eq[band].reset_coeffs_ext<fwd_idx> (
          co_biquad, freq, q, sr, highpass_tag {});
      }
      else {
        _eq[band].reset_coeffs_ext<fwd_idx> (
          co_biquad, freq, q, sr, biquad::mvic_highpass_hq_tag {});
      }
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
      reset_latency();
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
    double samplerate,
    double snr_db)
  {
    std::array<double_x1, biquad::n_coeffs> co;
    biquad::reset_coeffs<double_x1> (
      co, vec_set<1> (freq), vec_set<1> (q), samplerate, bandpass_tag {});
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

      freq[b] = get_parameter (freq_param {}).min;
      freq[b] = midi_note_to_hz (freq[b]);
      q[b]    = get_parameter (q_param {}).max;
    });

    for (uint b = 0; b < n_bands; ++b) {
      // Tuning was 11 stages on normal quality for the two first bands and 10
      // stages for the other two, both at 48000KHz. The lowest quality has a
      // truncation noise SNR of 90dB. The medium 180dB and the higher 360dB.
      _n_stages[b][0] = get_n_stages (freq[b], q[b], samplerate, 90.);
      _n_stages[b][1] = _n_stages[b][0] + 1;
      _n_stages[b][2] = _n_stages[b][1] + 1;
      static_assert (n_quality_steps == 3, "Update this!");

      for (auto& stages : _n_stages[b]) {
        stages = std::min (stages, max_n_stages);
      }
    }
  }
  //----------------------------------------------------------------------------
  void clear_orders()
  {
    for (auto& elem : _order_serial) {
      elem = -1;
    }
    for (auto& elem : _order_parallel) {
      elem = -1;
    }
  }
  //----------------------------------------------------------------------------
  enum class bandtype {
    off,
    bell,
    lshelf,
    hshelf,
    lp,
    hp,
    count,
  };
  //----------------------------------------------------------------------------
  bool is_parallel (bandtype t)
  {
    return t >= bandtype::bell && t <= bandtype::hshelf;
  }
  //----------------------------------------------------------------------------
  struct bandconfig {
    double   slope_gain       = 0.;
    bandtype type             = bandtype::off;
    float    freq_note        = constexpr_midi_note_to_hz (440.f);
    float    q                = 0.7f;
    float    gain_db          = 0.f;
    float    diff             = 0.f;
    uint     latency          = 0;
    bool     has_changes      = false;
    bool     reset_band_state = false;
  };
  //----------------------------------------------------------------------------
  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  enum class topology { stereo, l, r, count };
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
      make_max_stages_t_rev<t_rev_pole_pair, max_n_stages>,
      biquad,
      biquad>,
    double_x2,
    true>;

  using coeff_array = std::array<eqs::value_type, eqs::n_coeffs>;

  alignas (eqs::value_type) std::array<coeff_array, n_bands> _target_coeffs;
  std::array<eqs, n_bands> _eq;

  std::array<delay_compensated_buffer<double_x2>, n_bands>            _bandcomp;
  delay_compensated_buffer<double_x2>                                 _incomp;
  std::vector<double_x2, overaligned_allocator<double_x2, sse_bytes>> _mem;

  std::array<bandconfig, n_bands>                        _cfg;
  std::array<int, n_bands>                               _order_serial;
  std::array<int, n_bands>                               _order_parallel;
  uint                                                   _sample;
  double                                                 _smooth_coeff;
  uint                                                   _quality;
  std::array<std::array<uint, n_quality_steps>, n_bands> _n_stages;
  uint                                                   _latency     = 0;
  plugin_context*                                        _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
