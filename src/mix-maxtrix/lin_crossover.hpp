#pragma once

#include "artv-common/dsp/own/classes/linphase_iir_crossover.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <uint Bands>
class mixmaxtrix_linphase_iir_crossover {
private:
  enum class paramtype { mode, frequency, diff, out, quality };

public:
  template <uint Band, paramtype Type>
  struct param {
    static constexpr uint band  = Band;
    static constexpr uint ptype = (uint) Type;
  };
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type     = dsp_types::crossover;
  static constexpr bus_types bus_type     = bus_types::stereo;
  static constexpr uint      n_inputs     = 1;
  static constexpr uint      n_bands      = Bands;
  static constexpr uint      n_outputs    = n_bands;
  static constexpr uint      n_crossovers = n_bands - 1;

  static_assert (Bands <= 4, "Implement more bands");
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _cfg                          = decltype (_cfg) {};
    _plugcontext                  = &pc;
    _t_spl                        = 1.f / pc.get_sample_rate();
    static constexpr float snr_db = 120;

    // Ultra Bad SNR
    _q_n_stages[0] = _crossv.get_n_stages (15000.f, _t_spl, snr_db - 60.);
    // Bad SNR
    _q_n_stages[1] = _crossv.get_n_stages (15000.f, _t_spl, snr_db - 30.);
    // Ultra treble
    _q_n_stages[2] = _crossv.get_n_stages (15000.f, _t_spl, snr_db);
    // Treble
    _q_n_stages[3] = _crossv.get_n_stages (5000.f, _t_spl, snr_db);
    // High mids
    _q_n_stages[4] = _crossv.get_n_stages (2500.f, _t_spl, snr_db);
    // Mids
    _q_n_stages[5] = _crossv.get_n_stages (1000.f, _t_spl, snr_db);
    // Low Mids
    _q_n_stages[6] = _crossv.get_n_stages (400.f, _t_spl, snr_db);
    // High Lows
    _q_n_stages[7] = _crossv.get_n_stages (180.f, _t_spl, snr_db);
    // Lows
    _q_n_stages[8] = _crossv.get_n_stages (100.f, _t_spl, snr_db);
    // Ultra Lows
    _q_n_stages[9] = _crossv.get_n_stages (50.f, _t_spl, snr_db);
    // High quality
    _q_n_stages[10] = _crossv.get_n_stages (20.f, _t_spl, snr_db);
    // Ultra quality
    _q_n_stages[11] = _crossv.get_n_stages (20.f, _t_spl, snr_db + 40.);
    // Insane
    _q_n_stages[12] = _crossv.get_n_stages (20.f, _t_spl, snr_db + 180.);

    for (auto& bcfg : _cfg) {
      bcfg.n_stages = _q_n_stages[0]; // temporary setting
    }

    reset_n_stages();
    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    _crossv.tick<T> (outs, ins, samples);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::mode>, int v)
  {
    static_assert (band < n_crossovers, "");
    _cfg[band].mode = v;
    uint latency    = _crossv.get_latency();
    update_crossv (band);
    uint new_latency = _crossv.get_latency();
    if (latency != new_latency) {
      _plugcontext->set_delay_compensation (new_latency);
    }
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::mode>)
  {
    return choice_param (
      1, make_cstr_array ("Off", "12dB/Oct", "24dB/Oct", "48dB/Oct"), 15);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::frequency>, float v)
  {
    static_assert (band < n_bands, "");
    _cfg[band].freq = midi_note_to_hz (v);
    update_crossv (band);
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::frequency>)
  {
    if constexpr (band == 0) {
      return frequency_parameter (20.0, 20000.0, 180.0);
    }
    else if constexpr (band == 1) {
      return frequency_parameter (20.0, 20000.0, 690.0);
    }
    else if constexpr (band == 2) {
      return frequency_parameter (20.0, 20000.0, 3400.0);
    }
    else {
      return frequency_parameter (20.0, 20000.0, 440.0);
    }
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::diff>, float v)
  {
    static_assert (band < n_bands, "");
    _cfg[band].diff = v * (0.01 * 0.2);
    update_crossv (band);
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::diff>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::out>, int v)
  {
    // dummy, to be used by processor
    static_assert (band < n_bands, "");
    _cfg[band].out = v + 1;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::out>)
  {
    constexpr auto str_array = make_cstr_array (
      "bus2", "bus3", "bus4", "bus5", "bus6", "bus7", "bus8");
    constexpr uint n_future_choices = 16;

    if constexpr (band == 0) {
      return choice_param (1, str_array, n_future_choices);
    }
    else if constexpr (band == 1) {
      return choice_param (2, str_array, n_future_choices);
    }
    else if constexpr (band == 2) {
      return choice_param (3, str_array, n_future_choices);
    }
    else {
      return choice_param (0, str_array, n_future_choices);
    }
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::quality>, int v)
  {
    assert (v >= 0 && v < _q_n_stages.size());
    uint stages = _q_n_stages[v];
    if (_cfg[band].n_stages == stages) {
      return;
    }
    _cfg[band].n_stages = stages;
    reset_n_stages();
  }

  static constexpr uint n_quality_steps = 13;

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::quality>)
  {
    auto str_array = make_cstr_array (
      "15kHz Broken",
      "15kHz Bad",
      "15kHz",
      "5kHz",
      "2.5kHz",
      "1kHz",
      "400Hz",
      "180Hz",
      "100Hz",
      "50Hz",
      "20Hz",
      "20Hz Ultra",
      "20Hz Insane");

    static_assert (str_array.size() == n_quality_steps);

    constexpr uint n_future_choices = 16;

    if constexpr (band == 0) {
      return choice_param (7, str_array, n_future_choices);
    }
    else if constexpr (band == 1) {
      return choice_param (6, str_array, n_future_choices);
    }
    else if constexpr (band == 2) {
      return choice_param (4, str_array, n_future_choices);
    }
    else {
      return choice_param (8, str_array, n_future_choices);
    }
  }
  //----------------------------------------------------------------------------
  using band1_mode_tag = param<0, paramtype::mode>;
  using band2_mode_tag = param<1, paramtype::mode>;
  using band3_mode_tag = param<2, paramtype::mode>;

  using band1_freq_tag = param<0, paramtype::frequency>;
  using band2_freq_tag = param<1, paramtype::frequency>;
  using band3_freq_tag = param<2, paramtype::frequency>;

  using band1_diff_tag = param<0, paramtype::diff>;
  using band2_diff_tag = param<1, paramtype::diff>;
  using band3_diff_tag = param<2, paramtype::diff>;

  using band1_out_tag = param<0, paramtype::out>;
  using band2_out_tag = param<1, paramtype::out>;
  using band3_out_tag = param<2, paramtype::out>;

  using band1_quality_tag = param<0, paramtype::quality>;
  using band2_quality_tag = param<1, paramtype::quality>;
  using band3_quality_tag = param<2, paramtype::quality>;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    band1_freq_tag,
    band1_diff_tag,
    band1_mode_tag,
    band1_out_tag,
    band2_freq_tag,
    band2_diff_tag,
    band2_mode_tag,
    band2_out_tag,
    band3_freq_tag,
    band3_diff_tag,
    band3_mode_tag,
    band3_out_tag,
    band1_quality_tag,
    band2_quality_tag,
    band3_quality_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void reset_n_stages()
  {
    std::array<uint, n_crossovers> stages;
    for (uint i = 0; i < n_crossovers; ++i) {
      stages[i] = _cfg[i].n_stages;
    }
    uint latency = _crossv.get_latency();
    _crossv.reset (_plugcontext->get_sample_rate(), stages);
    for (uint i = 0; i < n_crossovers; ++i) {
      update_crossv (i);
    }
    uint new_latency = _crossv.get_latency();
    if (latency != new_latency) {
      _plugcontext->set_delay_compensation (new_latency);
    }
  }
  //----------------------------------------------------------------------------
  void update_crossv (uint n)
  {
    float diff = _cfg[n].diff;
    float f    = _cfg[n].freq;
    _crossv.set_crossover_point (
      n,
      f + (f * diff),
      f - (f * diff),
      ((uint) (_cfg[n].mode != 0)) << _cfg[n].mode);
  }
  //----------------------------------------------------------------------------
  struct bandcfg {
    float freq     = 30.f;
    float diff     = 0.f;
    uint  mode     = 0;
    uint  out      = 0;
    uint  n_stages = 0;
  };
  //----------------------------------------------------------------------------
  std::array<bandcfg, n_crossovers> _cfg;
  linphase_iir_crossover<3>         _crossv;
  float                             _snr_db = 120.;
  float                             _t_spl;
  plugin_context*                   _plugcontext = nullptr;
  std::array<uint, n_quality_steps> _q_n_stages;
};
//------------------------------------------------------------------------------
} // namespace artv
