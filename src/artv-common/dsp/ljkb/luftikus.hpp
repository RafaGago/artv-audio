#pragma once

#include <disthro-ports/ports-legacy/luftikus/source/dsp/eqdsp.h>

#include "artv-common/dsp/own/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace ljkb {
//------------------------------------------------------------------------------
class luftikus {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::eq;
  //----------------------------------------------------------------------------
  struct gain_10hz_tag {};
  void set (gain_10hz_tag, float v)
  {
    for (auto& eq : _eqs) {
      eq.setGain (eq_dsp::kBand10, v);
    }
  }
  static constexpr auto get_parameter (gain_10hz_tag)
  {
    return float_param ("dB", -10.0, 10.0, 0.0, 0.5);
  }
  //----------------------------------------------------------------------------
  struct gain_40hz_tag {};
  void set (gain_40hz_tag, float v)
  {
    for (auto& eq : _eqs) {
      eq.setGain (eq_dsp::kBand40, v);
    }
  }
  static constexpr auto get_parameter (gain_40hz_tag)
  {
    return float_param ("dB", -10.0, 10.0, 0.0, 0.5);
  }
  //----------------------------------------------------------------------------
  struct gain_160hz_tag {};
  void set (gain_160hz_tag, float v)
  {
    for (auto& eq : _eqs) {
      eq.setGain (eq_dsp::kBand160, v);
    }
  }
  static constexpr auto get_parameter (gain_160hz_tag)
  {
    return float_param ("dB", -10.0, 10.0, 0.0, 0.5);
  }
  //----------------------------------------------------------------------------
  struct gain_640hz_tag {};
  void set (gain_640hz_tag, float v)
  {
    for (auto& eq : _eqs) {
      eq.setGain (eq_dsp::kBand640, v);
    }
  }
  static constexpr auto get_parameter (gain_640hz_tag)
  {
    return float_param ("dB", -10.0, 10.0, 0.0, 0.5);
  }
  //----------------------------------------------------------------------------
  struct gain_2k5_tag {};
  void set (gain_2k5_tag, float v)
  {
    for (auto& eq : _eqs) {
      eq.setGain (eq_dsp::kShelf2k5, v);
    }
  }
  static constexpr auto get_parameter (gain_2k5_tag)
  {
    return float_param ("dB", -10.0, 10.0, 0.0, 0.5);
  }
  //----------------------------------------------------------------------------
  struct gain_hi_tag {};
  void set (gain_hi_tag, float v)
  {
    for (auto& eq : _eqs) {
      eq.setGain (eq_dsp::kShelfHi, v);
    }
  }
  static constexpr auto get_parameter (gain_hi_tag)
  {
    return float_param ("dB", -10.0, 10.0, 0.0, 0.5);
  }
  //----------------------------------------------------------------------------
  struct freq_hi_tag {};
  void set (freq_hi_tag, int v)
  {
    for (auto& eq : _eqs) {
      eq.setHighShelf ((eq_dsp::HighShelf) v);
    }
  }
  static constexpr auto get_parameter (freq_hi_tag)
  {
    return choice_param (
      0, make_cstr_array ("Off", "2.5kHz", "5kHz", "10kHz", "20kHz", "40kHz"));
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void set (mode_tag, int v)
  {
    for (auto& eq : _eqs) {
      eq.setMastering (!!(v & 1));
      eq.setAnalog (!!((v >> 1) & 1));
    }
  }
  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0, make_cstr_array ("Normal", "Mastering", "Analog", "Analog Mastering"));
  }
  //----------------------------------------------------------------------------
  struct keep_gain_tag {};
  void set (keep_gain_tag, int v)
  {
    for (auto& eq : _eqs) {
      eq.setKeepGain (!!v);
    }
  }
  static constexpr auto get_parameter (keep_gain_tag)
  {
    return choice_param (0, make_cstr_array ("Off", "On"));
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    gain_10hz_tag,
    gain_40hz_tag,
    gain_160hz_tag,
    gain_640hz_tag,
    gain_2k5_tag,
    gain_hi_tag,
    freq_hi_tag,
    mode_tag,
    keep_gain_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    for (auto& eq : _eqs) {
      eq.setBlockSize (pc.get_max_block_samples());
      eq.setSampleRate ((float) pc.get_sample_rate());
    }

    mp11::mp_for_each<parameters> ([=] (auto param) {
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  void process_block_replacing (std::array<float*, 2> chnls, uint block_samples)
  {
    _eqs[0].processBlock (chnls[0], block_samples);
    _eqs[1].processBlock (chnls[1], block_samples);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  using eq_dsp = ::artv_dsp_pull::luftikus::EqDsp;
  std::array<eq_dsp, 2> _eqs;
};
//------------------------------------------------------------------------------
}} // namespace artv::ljkb
