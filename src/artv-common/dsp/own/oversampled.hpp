#pragma once

#include "artv-common/dsp/own/delay_line.hpp"
#include "artv-common/dsp/own/fir.hpp"
#include "artv-common/dsp/own/misc.hpp"
#include "artv-common/dsp/own/oversampled_coeffs.hpp"
#include "artv-common/dsp/own/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
class oversampled_plugin_context : public plugin_context {
public:
  virtual ~oversampled_plugin_context() {};

  plugin_context*            pc                 = nullptr;
  uint                       oversampling_order = 0;
  std::function<uint (uint)> on_set_delay_compensation;

  uint get_oversampling() const { return 1u << oversampling_order; }

  uint get_sample_rate() const override
  {
    return pc->get_sample_rate() * get_oversampling();
  };

  uint get_max_block_samples() const override
  {
    return pc->get_max_block_samples();
  };

  plugin_play_state get_play_state() const override
  {
    return pc->get_play_state();
  }

  void set_delay_compensation (uint samples) override
  {
    pc->set_delay_compensation (on_set_delay_compensation (samples));
  };
};
//------------------------------------------------------------------------------
// an external adaptor class to oversample a DSP module.
static constexpr uint oversampled_max_oversampling = 16;

template <class fx, class T = float>
class oversampled {
public:
  static constexpr uint      max_oversampling = oversampled_max_oversampling;
  static constexpr uint      num_channels     = 2;
  static constexpr dsp_types dsp_type         = fx::dsp_type;
  using fx_type                               = fx;
  //----------------------------------------------------------------------------
  struct oversampling_tag {};

  void set (oversampling_tag, int v)
  {
    // design flaw, changing the oversample rate requires all the parameters to
    // be set manually. As on JUCE they aren't event based, but set once per
    // block a bad block after changing the sample rate is assumed.
    if (_pc.oversampling_order != v) {
      uint delay_prev        = get_total_delay();
      _pc.oversampling_order = v;
      _fx.reset (_pc);
      reset_predelay();
      reset_resamplers();
      uint delay = get_total_delay();
      if (delay_prev != delay) {
        _pc.pc->set_delay_compensation (delay);
      }
    }
  }

  static constexpr auto get_parameter (oversampling_tag)
  {
    return choice_param (
      0, make_cstr_array ("1x", "2x", "4x", "8x", "16x"), 30);
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp11::mp_push_back<typename fx::parameters, oversampling_tag>;
  //----------------------------------------------------------------------------
  template <class Tag, class U>
  void set (Tag t, U v)
  {
    _fx.set (t, v);
  }

  template <class Tag>
  static constexpr auto get_parameter (Tag t)
  {
    fx::get_parameter (t);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _pc.pc = &pc;
    _pc.on_set_delay_compensation
      = [=] (uint samples) { return delay_compensation_changing (samples); };
    _pc.oversampling_order = 0;
    _predelay_samples      = 0;
    _module_delay          = 0;
    _oversampler_delay     = 0;
    _predelay_line.reset (num_channels, max_oversampling);
    _fx.reset (_pc);
  }
  //----------------------------------------------------------------------------
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    assert (_work_buffer.size() >= (2 * samples * _pc.get_oversampling()));
    if (_pc.oversampling_order == 0) {
      // fast path
      _fx.process_block_replacing (chnls, samples);
      return;
    }
    // tell the compiler that these won't change during the loop.
    uint              predelay_samples   = _predelay_samples;
    uint              oversampling_order = _pc.oversampling_order;
    uint              ratio              = 1u << oversampling_order;
    std::array<T*, 2> upsampled
      = {&_work_buffer[0], &_work_buffer[_work_buffer.size() / 2]};

    for (uint i = 0; i < samples; ++i) {
      std::array<T, 2> interleaved {chnls[0][i], chnls[1][i]};
      if (predelay_samples) {
        _predelay_line.push (interleaved);
        _predelay_line.get (interleaved, predelay_samples);
      }
      uint idx = (i << oversampling_order);
      auto l   = make_crange (upsampled[0] + idx, ratio);
      auto r   = make_crange (upsampled[1] + idx, ratio);
      _interpolator.tick ({l, r}, interleaved);
    }

    _fx.process_block_replacing (upsampled, (samples << oversampling_order));

    for (uint i = 0; i < samples; ++i) {
      uint idx    = (i << oversampling_order);
      auto l      = make_crange (upsampled[0] + idx, ratio);
      auto r      = make_crange (upsampled[1] + idx, ratio);
      auto ret    = _decimator.tick ({l, r});
      chnls[0][i] = ret[0];
      chnls[1][i] = ret[1];
    }
  }
  //----------------------------------------------------------------------------
  // for cache friendliness this buffer is external.
  void set_oversampling_work_buffer (crange<T> buff) { _work_buffer = buff; }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  uint delay_compensation_changing (uint samples)
  {
    _module_delay = samples;
    reset_predelay();
    return get_total_delay();
  }
  //----------------------------------------------------------------------------
  template <uint ratio>
  void reset_oversamplers()
  {
    _oversampler_delay = linear_phase_fir_coeffs<ratio>::latency();
    _oversampler_delay *= 2; // Up and downsampler, 2 times the latency
    _interpolator.reset (linear_phase_fir_coeffs<ratio>::data(), ratio, true);
    _decimator.reset (linear_phase_fir_coeffs<ratio>::data(), ratio);
  }
  //----------------------------------------------------------------------------
  void reset_predelay()
  {
    _predelay_samples   = 0;
    uint next_int_delay = round_ceil (_module_delay, _pc.get_oversampling());
    // making the resulting latency an integer by adding util we reach a full
    // sample at the base sample rate
    _predelay_samples = next_int_delay - _module_delay;
  }
  //----------------------------------------------------------------------------
  void reset_resamplers()
  {
    switch (_pc.oversampling_order) {
    case 0:
      _oversampler_delay = 0;
      break;
    case 1:
      // TODO: halfband optimization
      reset_oversamplers<2>();
      break;
    case 2:
      reset_oversamplers<4>();
      break;
    case 3:
      reset_oversamplers<8>();
      break;
    case 4:
      reset_oversamplers<16>();
      break;
    default:
      break;
    }
  }
  //----------------------------------------------------------------------------
  uint get_total_delay() const
  {
    uint module_delay
      = (_module_delay + _predelay_samples) / _pc.get_oversampling();
    return _oversampler_delay + module_delay;
  }
  //----------------------------------------------------------------------------
  fx                           _fx;
  oversampled_plugin_context   _pc;
  delay_line<T>                _predelay_line;
  uint                         _predelay_samples;
  uint                         _module_delay;
  uint                         _oversampler_delay;
  uint                         _informed_delay;
  lth_band_fir_decimator<T, 2> _decimator;
  fir_interpolator<T, 2>       _interpolator;
  crange<T>                    _work_buffer;
};
//------------------------------------------------------------------------------

} // namespace artv
