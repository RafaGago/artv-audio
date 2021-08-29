#pragma once

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/fir.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/oversampled_coeffs.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/delay_compensation_buffers.hpp"
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
  static constexpr uint n_channels       = 2;

  using decimator_type    = lth_band_fir_decimator<T, 2, 16>;
  using interpolator_type = fir_interpolator<T, 2, 16>;
  //----------------------------------------------------------------------------
  void set_oversampling_order (
    uint                                  order,
    std::function<void (plugin_context&)> fx_reset_fn,
    interpolator_type&                    interpolator,
    decimator_type*                       decimator)
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
    _n_samples_fractional  = 0;
    _module_delay          = 0;
    _oversample_delay      = 0;
    // the +1 there is to do fractional delay compensations.
    uint bfsz = (pc.get_max_block_samples() + 1) * max_oversampling;
    _work_buffers_mem.clear();
    _work_buffers_mem.resize (bfsz * n_channels);
    for (uint i = 0; i < n_channels; ++i) {
      _work_buffers[i].reset (
        {&_work_buffers_mem[i * bfsz], bfsz}, max_oversampling);
    }
    fx_reset_fn (_pc);
  }
  //----------------------------------------------------------------------------
  void upsample (std::array<T*, 2> chnls, uint samples, interpolator_type& up)
  {}
  //----------------------------------------------------------------------------
  template <bool downsample>
  void process_block_replacing (
    std::array<T*, 2>                             chnls,
    uint                                          samples,
    std::function<void (std::array<T*, 2>, uint)> fx_process_block_replacing,
    interpolator_type&                            interpolator,
    decimator_type*                               decimator)
  {
    if (_pc.oversampling_order == 0) {
      // fast path, oversampling disabled
      fx_process_block_replacing (chnls, samples);
      return;
    }

    // tell the compiler that these won't change during the loop.
    uint n_frac_spl         = _n_samples_fractional;
    uint oversampling_order = _pc.oversampling_order;
    uint ratio              = 1u << oversampling_order;

    // the head of the work buffer is reserved for storing previous samples for
    // fractional delay corrections.
    std::array<T*, 2> upsampled = {
      _work_buffers[0].get_write_buffer().data(),
      _work_buffers[1].get_write_buffer().data()};

    for (uint i = 0; i < samples; ++i) {
      std::array<T, 2> interleaved {chnls[0][i], chnls[1][i]};
      uint             idx = (i << oversampling_order);
      auto             l   = make_crange (upsampled[0] + idx, ratio);
      auto             r   = make_crange (upsampled[1] + idx, ratio);
      interpolator.tick ({l, r}, interleaved);
    }
    // The fractional delay line is placed after the upsampler
    std::array<T*, 2> upsampled_frac_corrected = {
      _work_buffers[0].get_read_buffer (n_frac_spl).data(),
      _work_buffers[1].get_read_buffer (n_frac_spl).data()};

    fx_process_block_replacing (
      upsampled_frac_corrected, (samples << oversampling_order));

    for (uint i = 0; i < samples; ++i) {
      uint idx = (i << oversampling_order);
      if constexpr (downsample) {
        auto l      = make_crange (upsampled_frac_corrected[0] + idx, ratio);
        auto r      = make_crange (upsampled_frac_corrected[1] + idx, ratio);
        auto ret    = decimator->tick ({l, r});
        chnls[0][i] = ret[0];
        chnls[1][i] = ret[1];
      }
      else {
        // just drop samples, on this configuration "fx_process_block_replacing"
        // is assumed to not introduce aliasing
        chnls[0][i] = upsampled_frac_corrected[0][idx];
        chnls[1][i] = upsampled_frac_corrected[1][idx];
      }
    }
    for (uint c = 0; c < n_channels; ++c) {
      _work_buffers[c].iteration_end (
        n_frac_spl, samples << oversampling_order);
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
  void reset_oversamplers (interpolator_type& up, decimator_type* down)
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
    _n_samples_fractional = 0;
    uint next_int_delay   = round_ceil (_module_delay, _pc.get_oversampling());
    // making the resulting latency an integer by adding util we reach a full
    // sample at the base sample rate
    _n_samples_fractional = next_int_delay - _module_delay;
  }
  //----------------------------------------------------------------------------
  void reset_resamplers (interpolator_type& up, decimator_type* down)
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
      = (_module_delay + _n_samples_fractional) / _pc.get_oversampling();
    return _oversample_delay + module_delay;
  }
  //----------------------------------------------------------------------------
  oversampled_plugin_context _pc;
  uint                       _n_samples_fractional;
  uint                       _module_delay;
  uint                       _oversample_delay;
  std::array<owned_delay_compensation_buffer<T>, n_channels> _work_buffers;
  std::vector<T>                                             _work_buffers_mem;
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
private:
  //----------------------------------------------------------------------------
  using interpolator_type = typename oversampled_common<T>::interpolator_type;
  using decimator_type    = typename oversampled_common<T>::decimator_type;

  using maybe_decimator
    = std::conditional_t<downsample, decimator_type, std::monostate>;

  fx                    _fx;
  oversampled_common<T> _common;
  interpolator_type     _interpolator;
  maybe_decimator       _decimator;
};
//------------------------------------------------------------------------------
template <class fx, class T = float>
using updownsampled = oversampled<fx, std::true_type, T>;

template <class fx, class T = float>
using upsampled = oversampled<fx, std::false_type, T>;

} // namespace artv
