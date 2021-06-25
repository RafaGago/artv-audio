#pragma once
// An FX that just upsamples and downsamples. It was just useful when developing

#include "artv-common/dsp/own/polyphase-fir.hpp"

namespace artv {

static std::array<double, 65> lowpass_coeffs = {
  // linphase_fir_sr88200hz_fc20000hz_tw4250hz_att96db = {
  5.057877639028037e-06,   2.9039377894283255e-06,  -3.3164365183698734e-05,
  -3.0727657843876072e-05, 9.4352090925605648e-05,  0.00013162273866776532,
  -0.00017880572133587201, -0.00038215090945201873, 0.00022485859781025392,
  0.00086973917326166718,  -8.7232733578806621e-05, -0.0016487139384029524,
  -0.00048015878750228816, 0.0026715301207921849,   0.001808844991548791,
  -0.0037084255323752797,  -0.0042562500115399787,  0.0042771152228502224,
  0.0081051094232479862,   -0.0036012320006442184,  -0.013443834577344457,
  0.00059225125635217789,  0.020065140542363375,    0.0062146260857539977,
  -0.027424172382569455,   -0.019021523565952866,   0.034685449744363023,
  0.042414161911683144,    -0.040862034229607332,   -0.092305872313849235,
  0.045018691964856743,    0.3135249364961305,      0.45351581310229289,
  0.3135249364961305,      0.045018691964856743,    -0.092305872313849235,
  -0.040862034229607332,   0.042414161911683144,    0.034685449744363023,
  -0.019021523565952866,   -0.027424172382569455,   0.0062146260857539977,
  0.020065140542363375,    0.00059225125635217789,  -0.013443834577344457,
  -0.0036012320006442184,  0.0081051094232479862,   0.0042771152228502224,
  -0.0042562500115399787,  -0.0037084255323752797,  0.001808844991548791,
  0.0026715301207921849,   -0.00048015878750228816, -0.0016487139384029524,
  -8.7232733578806621e-05, 0.00086973917326166718,  0.00022485859781025392,
  -0.00038215090945201873, -0.00017880572133587201, 0.00013162273866776532,
  9.4352090925605648e-05,  -3.0727657843876072e-05, -3.3164365183698734e-05,
  2.9039377894283255e-06,  5.057877639028037e-06,
};
//------------------------------------------------------------------------------
class polyphase_fir_test {
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
    _up.reset (lowpass_coeffs, 2, 2);
    _down.reset (lowpass_coeffs, 2, 2);
    _gain = 1.f;
  }
  //----------------------------------------------------------------------------
  void process_block_replacing (std::array<float*, 2> chnls, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {
      std::array<std::array<double, 2>, 2> upsampled;
      _up.tick (upsampled[0], chnls[0][i], 0);
      _up.tick (upsampled[1], chnls[1][i], 1);
      chnls[0][i] = _down.tick (upsampled[0], 0) * _gain;
      chnls[1][i] = _down.tick (upsampled[1], 1) * _gain;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  fir_decimator<double>    _down;
  fir_interpolator<double> _up;
  float                    _gain = 1.;
  //----------------------------------------------------------------------------
};

} // namespace artv
