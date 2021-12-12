#pragma once

#include "artv-common/dsp/own/classes/linphase_iir_crossover.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <uint Bands>
class mixmaxtrix_linphase_iir_crossover {
private:
  enum class paramtype { mode, frequency, diff, out };

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
    _cfg         = decltype (_cfg) {};
    _plugcontext = &pc;
    _snr_db      = 120.;
    reset_snr();

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    _plugcontext->set_delay_compensation (_crossv.get_latency());
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
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
      // TODO: latency change. Will this be propagated for crossovers?
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
      return frequency_parameter (50.0, 20000.0, 140.0);
    }
    else if constexpr (band == 1) {
      return frequency_parameter (120.0, 20000.0, 390.0);
    }
    else if constexpr (band == 2) {
      return frequency_parameter (250.0, 20000.0, 2200.0);
    }
    else {
      return frequency_parameter (400.0, 20000.0, 440.0);
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
  struct snr_tag {};
  void set (snr_tag, int v)
  {
    uint snr = v * 20;
    snr += 80;
    if ((float) snr == _snr_db) {
      return;
    }
    _snr_db = snr;
    reset_snr();
    for (uint i = 0; i < n_crossovers; ++i) {
      update_crossv (i);
    }
    _plugcontext->set_delay_compensation (_crossv.get_latency());
  }

  static constexpr auto get_parameter (snr_tag)
  {
    return choice_param (
      2,
      make_cstr_array (
        "80dB",
        "100dB",
        "120dB",
        "140dB",
        "160dB",
        "180dB",
        "200dB",
        "220dB",
        "240dB",
        "260dB"),
      16);
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
    snr_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void reset_snr()
  {
    std::array<float, n_crossovers> min_freq;
    mp11::mp_for_each<mp11::mp_iota_c<n_crossovers>> ([&] (auto i) {
      static constexpr uint idx = decltype (i)::value;
      auto note   = get_parameter (param<idx, paramtype::frequency> {}).min;
      min_freq[i] = midi_note_to_hz (note);
    });
    _crossv.reset (_plugcontext->get_sample_rate(), min_freq, _snr_db);
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
    float freq = 30.f;
    float diff = 0.f;
    uint  mode = 0;
    uint  out  = 0;
  };
  //----------------------------------------------------------------------------
  std::array<bandcfg, n_crossovers> _cfg;
  linphase_iir_crossover<3>         _crossv;
  float                             _snr_db      = 120.;
  plugin_context*                   _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv