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
// Utility class implementing oversampling methods. Its only purpose is to allow
// the compiler to reduce code bloat by not implementing everything on the
// "oversampled" class below, I guess that LTO and compile flags could get rid
// of this, doing so for not having to remember to check.
template <class T = float>
class oversampled_common {
public:
  static constexpr uint max_oversampling = 16;
  static constexpr uint num_channels     = 2;
  //----------------------------------------------------------------------------
  void set_oversampling_order (
    uint                                  order,
    std::function<void (plugin_context&)> fx_reset_fn,
    fir_interpolator<T, 2>&               interpolator,
    lth_band_fir_decimator<T, 2>*         decimator)
  {
    if (_pc.oversampling_order != order) {
      uint delay_prev        = get_total_delay();
      _pc.oversampling_order = order;
      fx_reset_fn (_pc);
      reset_predelay();
      reset_resamplers (interpolator, decimator);
      uint delay = get_total_delay();
      if (delay_prev != delay) {
        _pc.pc->set_delay_compensation (delay);
      }
    }
  }
  //----------------------------------------------------------------------------
  void reset (
    plugin_context&                       pc,
    std::function<void (plugin_context&)> fx_reset_fn)
  {
    _pc.pc = &pc;
    _pc.on_set_delay_compensation
      = [=] (uint samples) { return delay_compensation_changing (samples); };
    _pc.oversampling_order = 0;
    _predelay_samples      = 0;
    _module_delay          = 0;
    _oversample_delay      = 0;
    _predelay_line.reset (num_channels, max_oversampling);
    fx_reset_fn (_pc);
  }
  //----------------------------------------------------------------------------
  void set_oversampling_work_buffer (crange<T> buff) { _work_buffer = buff; }
  //----------------------------------------------------------------------------
  void upsample (
    std::array<T*, 2>       chnls,
    uint                    samples,
    fir_interpolator<T, 2>& up)
  {}
  //----------------------------------------------------------------------------
  template <bool downsample>
  void process_block_replacing (
    std::array<T*, 2>                             chnls,
    uint                                          samples,
    std::function<void (std::array<T*, 2>, uint)> fx_process_block_replacing,
    fir_interpolator<T, 2>&                       interpolator,
    lth_band_fir_decimator<T, 2>*                 decimator)
  {
    assert (_work_buffer.size() >= (2 * samples * _pc.get_oversampling()));
    if (_pc.oversampling_order == 0) {
      // fast path
      fx_process_block_replacing (chnls, samples);
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
      interpolator.tick ({l, r}, interleaved);
    }

    fx_process_block_replacing (upsampled, (samples << oversampling_order));

    for (uint i = 0; i < samples; ++i) {
      uint idx = (i << oversampling_order);
      if constexpr (downsample) {
        auto l      = make_crange (upsampled[0] + idx, ratio);
        auto r      = make_crange (upsampled[1] + idx, ratio);
        auto ret    = decimator->tick ({l, r});
        chnls[0][i] = ret[0];
        chnls[1][i] = ret[1];
      }
      else {
        // just drop samples, the process is assumed to not introduce aliasing
        chnls[0][i] = upsampled[0][idx];
        chnls[1][i] = upsampled[1][idx];
      }
    }
  }
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
  void reset_oversamplers (
    fir_interpolator<T, 2>&       up,
    lth_band_fir_decimator<T, 2>* down)
  {
    _oversample_delay = linear_phase_fir_coeffs<ratio>::latency();
    up.reset (linear_phase_fir_coeffs<ratio>::data(), ratio, true);
    if (down) {
      _oversample_delay *= 2; // Up and downsampler, 2 times the latency
      down->reset (linear_phase_fir_coeffs<ratio>::data(), ratio);
    }
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
  void reset_resamplers (
    fir_interpolator<T, 2>&       up,
    lth_band_fir_decimator<T, 2>* down)
  {
    switch (_pc.oversampling_order) {
    case 0:
      _oversample_delay = 0;
      break;
    case 1:
      // TODO: halfband optimization
      reset_oversamplers<2> (up, down);
      break;
    case 2:
      reset_oversamplers<4> (up, down);
      break;
    case 3:
      reset_oversamplers<8> (up, down);
      break;
    case 4:
      reset_oversamplers<16> (up, down);
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
    return _oversample_delay + module_delay;
  }
  //----------------------------------------------------------------------------
  oversampled_plugin_context _pc;
  delay_line<T>              _predelay_line;
  uint                       _predelay_samples;
  uint                       _module_delay;
  uint                       _oversample_delay;
  crange<T>                  _work_buffer;
};
//------------------------------------------------------------------------------
// an external adaptor class to oversample a DSP module.
// notice that "mode" gets a tag class so "is_same_template" can be used
// (reminder: non-typename template parameters break variadic template args).
struct oversampled_amount_tag {};

template <class fx, class downsamples = std::true_type, class T = float>
class oversampled {
public:
  static constexpr bool      downsample = downsamples::value;
  static constexpr dsp_types dsp_type   = fx::dsp_type;
  using fx_type                         = fx;
  //----------------------------------------------------------------------------
  void set (oversampled_amount_tag, int v)
  {
    // design flaw(ish), changing the sample rate requires a module reset, which
    // leaves the parameters on its default value. The modules have getters but
    // not setters for the parameters. The controlling class has the unwritten
    // responsibility to reset all the parameters the effects after a sample
    // rate change when using this convenience class.
    auto fx_reset = [=] (plugin_context& pc) { _fx.reset (pc); };
    if constexpr (downsample) {
      _common.set_oversampling_order (v, fx_reset, _interpolator, &_decimator);
    }
    else {
      _common.set_oversampling_order (v, fx_reset, _interpolator, nullptr);
    }
  }

  static constexpr auto get_parameter (oversampled_amount_tag)
  {
    return choice_param (
      0, make_cstr_array ("1x", "2x", "4x", "8x", "16x"), 30);
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp11::mp_push_back<typename fx::parameters, oversampled_amount_tag>;
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
    auto fx_reset = [=] (plugin_context& pc) { _fx.reset (pc); };
    _common.reset (pc, fx_reset);
  }
  //----------------------------------------------------------------------------
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    auto fx_process = [=] (std::array<T*, 2> c, uint s) {
      _fx.process_block_replacing (c, s);
    };
    if constexpr (downsample) {
      _common.template process_block_replacing<downsample> (
        chnls, samples, fx_process, _interpolator, &_decimator);
    }
    else {
      _common.template process_block_replacing<downsample> (
        chnls, samples, fx_process, _interpolator, nullptr);
    }
  }
  //----------------------------------------------------------------------------
  // for cache friendliness this buffer is external.
  void set_oversampling_work_buffer (crange<T> buff)
  {
    _common.set_oversampling_work_buffer (buff);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  using maybe_decimator = std::
    conditional_t<downsample, lth_band_fir_decimator<T, 2>, std::monostate>;

  fx                     _fx;
  oversampled_common<T>  _common;
  fir_interpolator<T, 2> _interpolator;
  maybe_decimator        _decimator;
};
//------------------------------------------------------------------------------
template <class fx, class T = float>
using updownsampled = oversampled<fx, std::true_type, T>;

template <class fx, class T = float>
using upsampled = oversampled<fx, std::false_type, T>;

} // namespace artv
