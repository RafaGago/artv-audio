#pragma once
// An FX that just upsamples and downsamples. It was just useful when developing

#include "artv-common/dsp/own/classes/fir.hpp"
#include "artv-common/dsp/own/classes/oversampled_coeffs.hpp"

namespace artv {

//------------------------------------------------------------------------------
class polyphase_fir_test {
private:
  static constexpr uint ratio    = 2;
  static constexpr uint channels = 2;
  using sample_type              = float;

public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::other;
  //----------------------------------------------------------------------------
  struct gain_tag {};

  void set (gain_tag, float v) { _gain = db_to_gain (v); }

  static constexpr auto get_parameter (gain_tag)
  {
    return float_param ("dB", -15.0, 15., 0.0, 0.25);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<gain_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _up.reset (linear_phase_fir_coeffs<ratio>::data(), ratio, true);
    _down.reset (linear_phase_fir_coeffs<ratio>::data(), ratio);
    _gain = 1.f;
  }
  //----------------------------------------------------------------------------
  void process_block_replacing (std::array<float*, 2> chnls, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {
      std::array<sample_type, channels> in = {chnls[0][i], chnls[1][i]};
      std::array<std::array<sample_type, ratio>, channels> upsampled;
      _up.tick ({make_crange (upsampled[0]), make_crange (upsampled[1])}, in);
#if 1
      auto ret
        = _down.tick ({make_crange (upsampled[0]), make_crange (upsampled[1])});
      chnls[0][i] = ret[0] * _gain;
      chnls[1][i] = ret[1] * _gain;
#else
      // only test upsampler...
      chnls[0][i] = upsampled[0][0] * _gain;
      chnls[1][i] = upsampled[1][0] * _gain;
#endif
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  lth_band_fir_decimator<sample_type, channels> _down;
  fir_interpolator<sample_type, channels>       _up;
  float                                         _gain = 1.;
  //----------------------------------------------------------------------------
};

} // namespace artv
