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
  struct out_diffusion_tag {};
  void set (out_diffusion_tag, float v)
  {
    if (v == _out_diffusion) {
      return;
    }
    _out_diffusion = v;
    _impl.set_output_diffusor_gain (v * 0.01f);
  }

  static constexpr auto get_parameter (out_diffusion_tag)
  {
    return float_param ("%", 0.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct early_gain_tag {};
  void set (early_gain_tag, float v)
  {
    if (v == _early_gain) {
      return;
    }
    _early_gain = v;
    _impl.set_early_gain (db_to_gain (v, -60.f));
  }

  static constexpr auto get_parameter (early_gain_tag)
  {
    return float_param ("dB", -60.f, 6.f, -6.f, 0.2f);
  }
  //----------------------------------------------------------------------------
  struct late_gain_tag {};
  void set (late_gain_tag, float v)
  {
    if (v == _late_gain) {
      return;
    }
    _late_gain = v;
    _impl.set_late_gain (db_to_gain (v, -60.f));
  }

  static constexpr auto get_parameter (late_gain_tag)
  {
    return float_param ("dB", -60.f, 6.f, -6.f, 0.2f);
  }
  //----------------------------------------------------------------------------
  struct early_2_late_bal_tag {};
  void set (early_2_late_bal_tag, float v)
  {
    if (v == _early_2_late_bal) {
      return;
    }
    _early_2_late_bal = v;
    v *= 0.01;
    _impl.set_in_to_late (1.f - v);
    _impl.set_early_to_late (v);
  }

  static constexpr auto get_parameter (early_2_late_bal_tag)
  {
    return float_param ("dB", 0.f, 100.f, -0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct time_tag {};
  void set (time_tag, float v)
  {
    if (v == _time_msec) {
      return;
    }
    _time_msec = v;
    _impl.set_time_msec (v);
  }

  static constexpr auto get_parameter (time_tag)
  {
    return float_param ("msec", 10.f, 20000.f, 1500.f, 1.f, 0.32f);
  }
  //----------------------------------------------------------------------------
  struct size_tag {};
  void set (size_tag, float v)
  {
    if (v == _size) {
      return;
    }
    _size = v;
    _impl.set_size (v * 0.01f);
  }

  static constexpr auto get_parameter (size_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.001f);
  }
  //----------------------------------------------------------------------------
  struct er_size_tag {};
  void set (er_size_tag, float v)
  {
    if (v == _er_size) {
      return;
    }
    _er_size = v;
    _impl.set_er_size (v * 0.01f);
  }

  static constexpr auto get_parameter (er_size_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.001f);
  }
  //----------------------------------------------------------------------------
  struct mod_freq_tag {};
  void set (mod_freq_tag, float v)
  {
    if (v == _mod_freq) {
      return;
    }
    _mod_freq = v;
    _impl.set_mod_freq (v * 0.01f);
  }

  static constexpr auto get_parameter (mod_freq_tag)
  {
    return float_param ("%", 0.f, 100.f, 5.f, 0.001f, 0.6f);
  }
  //----------------------------------------------------------------------------
  struct mod_depth_tag {};
  void set (mod_depth_tag, float v)
  {
    if (v == _mod_depth) {
      return;
    }
    _mod_depth = v;
    _impl.set_mod_depth (v * 0.01f);
  }

  static constexpr auto get_parameter (mod_depth_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct mod_stereo_tag {};
  void set (mod_stereo_tag, float v)
  {
    if (v == _mod_stereo) {
      return;
    }
    _mod_stereo = v;
    _impl.set_mod_stereo (v * 0.01f);
  }

  static constexpr auto get_parameter (mod_stereo_tag)
  {
    return float_param ("%", -100.f, 100.f, -50.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct mod_wave_tag {};
  void set (mod_wave_tag, uint v)
  {
    if (v == _mod_wave) {
      return;
    }
    _mod_wave = v;
    _impl.set_mod_wave (v);
  }

  static constexpr auto get_parameter (mod_wave_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Sample&Hold", "Sine", "Triangle", "Square", "Saw", "Trapezoid"),
      10);
  }
  //----------------------------------------------------------------------------
  struct l_matrix_angle_tag {};
  void set (l_matrix_angle_tag, float v)
  {
    if (v == _l_matrix_angle) {
      return;
    }
    _l_matrix_angle = v;
    _impl.set_l_matrix_angle (v * 0.01f);
  }

  static constexpr auto get_parameter (l_matrix_angle_tag)
  {
    return float_param ("%", 0.f, 100.f, 52.34f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct r_matrix_angle_tag {};
  void set (r_matrix_angle_tag, float v)
  {
    if (v == _r_matrix_angle) {
      return;
    }
    _r_matrix_angle = v;
    _impl.set_r_matrix_angle (v * 0.01f);
  }

  static constexpr auto get_parameter (r_matrix_angle_tag)
  {
    return float_param ("%", 0.f, 100.f, 48.12f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct lr_matrix_angle_tag {};
  void set (lr_matrix_angle_tag, float v)
  {
    if (v == _lr_matrix_angle) {
      return;
    }
    _lr_matrix_angle = v;
    _impl.set_lr_matrix_angle (v * 0.01f);
  }

  static constexpr auto get_parameter (lr_matrix_angle_tag)
  {
    return float_param ("%", 0.f, 100.f, 51.23f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct damp_freq_tag {};
  void set (damp_freq_tag, float v)
  {
    if (v == _damp_freq) {
      return;
    }
    _damp_freq = v;
    _impl.set_damp_freq (_damp_freq * 0.01f);
  }

  static constexpr auto get_parameter (damp_freq_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct damp_factor_tag {};
  void set (damp_factor_tag, float v)
  {
    if (v == _damp_factor) {
      return;
    }
    _damp_factor = v;
    _impl.set_damp_factor (_damp_factor * 0.01f);
  }

  static constexpr auto get_parameter (damp_factor_tag)
  {
    return float_param ("%", 0.f, 100.f, 40.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct hp_freq_tag {};
  void set (hp_freq_tag, float v)
  {
    if (v == _hp_freq) {
      return;
    }
    _hp_freq = v;
    _impl.set_hp_freq (_hp_freq * 0.01f);
  }

  static constexpr auto get_parameter (hp_freq_tag)
  {
    return float_param ("%", 0.f, 100.f, 35.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct lf_time_factor_tag {};
  void set (lf_time_factor_tag, float v)
  {
    if (v == _lf_rt60_factor) {
      return;
    }
    _lf_rt60_factor = v;
    _impl.set_lf_time_factor (_lf_rt60_factor * 0.01f);
  }

  static constexpr auto get_parameter (lf_time_factor_tag)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct predelay_tag {};
  void set (predelay_tag, float v)
  {
    if (v == _predelay) {
      return;
    }
    _predelay = v;
    _impl.set_predelay (v);
  }

  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("sixteenths", 0.f, 16.f, 0.f, 0.001f);
  }
  //----------------------------------------------------------------------------
  struct gap_tag {};
  void set (gap_tag, float v)
  {
    if (v == _gap) {
      return;
    }
    _gap = v;
    _impl.set_gap (v);
  }

  static constexpr auto get_parameter (gap_tag)
  {
    return float_param ("sixteenths", 0.f, 16.f, 0.f, 0.001f);
  }
#if 0
  //----------------------------------------------------------------------------
  struct test_param_tag {};
  void set (test_param_tag, uint v) { _impl.set_test_param (v); }

  static constexpr auto get_parameter (test_param_tag)
  {
    return int_param ("test", 0, 10, 0);
  }
#endif
  //----------------------------------------------------------------------------
#if 0
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
    _impl.reset (
      pc.get_sample_rate(),
      pc.get_play_state().bpm,
      _impl.get_default_cfg_preset());

    _impl.set_mod_stereo (-0.5f);

    static constexpr float invalid_val = INFINITY;
    _time_msec                         = invalid_val;
    _size                              = invalid_val;
    _er_size                           = invalid_val;
    _in_diffusion                      = invalid_val;
    _out_diffusion                     = invalid_val;
    _early_gain                        = invalid_val;
    _late_gain                         = invalid_val;
    _early_2_late_bal                  = invalid_val;
    _mod_freq                          = invalid_val;
    _mod_depth                         = invalid_val;
    _mod_stereo                        = invalid_val;
    _l_matrix_angle                    = invalid_val;
    _r_matrix_angle                    = invalid_val;
    _lr_matrix_angle                   = invalid_val;
    _damp_freq                         = invalid_val;
    _damp_factor                       = invalid_val;
    _hp_freq                           = invalid_val;
    _lf_rt60_factor                    = invalid_val;
    _predelay                          = invalid_val;
    _gap                               = invalid_val;
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

  float _time_msec;
  float _size;
  float _er_size;
  float _in_diffusion;
  float _out_diffusion;
  float _early_gain;
  float _late_gain;
  float _early_2_late_bal;
  float _mod_freq;
  float _mod_depth;
  float _mod_stereo;
  float _l_matrix_angle;
  float _r_matrix_angle;
  float _lr_matrix_angle;
  float _damp_freq;
  float _damp_factor;
  float _hp_freq;
  float _lf_rt60_factor;
  float _predelay;
  float _gap;
  uint  _mod_wave;
};
//------------------------------------------------------------------------------
} // namespace artv
