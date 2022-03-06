#pragma once

#pragma once

#include <numeric>

#include "artv-common/dsp/own/classes/fdn_stereo_8.hpp"

namespace artv {

//------------------------------------------------------------------------------
class reverb {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct in_diffusion_tag {};
  void set (in_diffusion_tag, float v)
  {
    if (v == _in_diffusion) {
      return;
    }
    _in_diffusion = v;
    _impl.set_input_diffusor_gain (v * 0.01f);
  }

  static constexpr auto get_parameter (in_diffusion_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct _early_gain {};
  void set (_early_gain, float v)
  {
    if (v == _early_att) {
      return;
    }
    _early_att = v;
    _impl.set_early_gain (db_to_gain (v, -60.f));
  }

  static constexpr auto get_parameter (_early_gain)
  {
    return float_param ("dB", -60.f, 6.f, -6.f, 0.2f);
  }
#if 0
  struct time_tag {};
  void set (time_tag, float v)
  {
    if (v == _fb_rt60_sec) {
      return;
    }
    _fb_rt60_sec = v;
  }

  static constexpr auto get_parameter (time_tag)
  {
    return float_param ("sec", 0.1f, 10.9f, 2.f, 0.01f, 0.5f);
  }
  //----------------------------------------------------------------------------
  struct predelay_tag {};
  void set (predelay_tag, float v) {}

  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("sec", 0.01f, 2.f, 0.1f, 0.01f, 0.6f);
  }
  //----------------------------------------------------------------------------
  struct er_late_tag {};
  void set (er_late_tag, float v) {}

  static constexpr auto get_parameter (er_late_tag)
  {
    return float_param ("%", 0.f, 100.f, 20.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct algorithm_tag {};
  void set (algorithm_tag, uint v)
  {
    if (v == _diff_mode) {
      return;
    }
    _diff_mode = v;
    reset_diffusor_times();
  }

  static constexpr auto get_parameter (algorithm_tag)
  {
    return choice_param (
      0, make_cstr_array ("1", "2", "3", "4", "5", "6", "7"), 64);
  }
  //----------------------------------------------------------------------------
  struct algorithm_param_tag {};
  void set (algorithm_param_tag, float v) {}

  static constexpr auto get_parameter (algorithm_param_tag)
  {
    return float_param ("%", 0.f, 100.f, 520.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct tail_gap_tag {};
  void set (tail_gap_tag, float v) {}

  static constexpr auto get_parameter (tail_gap_tag)
  {
    return float_param ("sec", 0.01f, 2.f, 0.1f, 0.01f, 0.6f);
  }
  //----------------------------------------------------------------------------
  struct tilt_tag {};
  void set (tilt_tag, float v) {}

  static constexpr auto get_parameter (tilt_tag)
  {
    return float_param ("%", 0.f, 100.f, 20.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct width_tag {};
  void set (width_tag, float v) {}

  static constexpr auto get_parameter (width_tag)
  {
    return float_param ("%", 0.f, 100.f, 20.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    time_tag,
    predelay_tag,
    er_late_tag,
    algorithm_tag,
    algorithm_param_tag,
    tail_gap_tag,
    tilt_tag,
    width_tag>;
#endif
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _impl.reset (pc.get_sample_rate(), _impl.get_default_cfg_preset());
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    _impl.process (outs, ins, samples);
  }
  //----------------------------------------------------------------------------
private:
  fdn_stereo_8 _impl;

  float _fb_rt60_sec;
  float _in_diffusion;
  float _early_att;
};
//------------------------------------------------------------------------------
} // namespace artv
