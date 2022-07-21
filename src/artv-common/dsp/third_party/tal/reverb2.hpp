#pragma once

#include <algorithm>
#include <array>
#include <optional>

#include <disthro-ports/ports-legacy/tal-reverb-2/source/Engine/ReverbEngine.h>

#include "artv-common/dsp/own/classes/add_ducker.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv { namespace tal {
//------------------------------------------------------------------------------
class reverb2 : private add_ducker<double_x2> {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct decay_time_tag {};
  void set (decay_time_tag, float v) { _reverb->setDecayTime (v * 0.01); }
  static constexpr auto get_parameter (decay_time_tag)
  {
    return float_param ("%", 0.0, 100.0, 20.0, 0.1);
  }
  //----------------------------------------------------------------------------
  struct predelay_tag {};
  void set (predelay_tag, float v)
  {
    _predelay = v;
    update_predelay();
  }
  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("ms", 0.0, 1000.0, 0.0, 0.01, 0.4);
  }
  //----------------------------------------------------------------------------
  struct predelay_sync_tag {};
  void set (predelay_sync_tag, float v)
  {
    _predelay_sync = v;
    update_predelay();
  }
  static constexpr auto get_parameter (predelay_sync_tag)
  {
    return float_param ("ms", 0.0, 1000.0, 0.0, 0.01, 0.4);
  }
  //----------------------------------------------------------------------------
  struct lowshelf_frequency_tag {};
  void set (lowshelf_frequency_tag, float v)
  {
    v = midi_note_to_hz (v);
    _reverb->setLowShelfFrequency (((v - 100.f) * (1.f / 9900.f)));
  }
  static constexpr auto get_parameter (lowshelf_frequency_tag)
  {
    return frequency_parameter (100.0, 10000.0, 200.0);
  }
  //----------------------------------------------------------------------------
  struct peak_frequency_tag {};
  void set (peak_frequency_tag, float v)
  {
    v = midi_note_to_hz (v);
    _reverb->setPeakFrequency ((v - 100.f) * (1.f / 9900.f));
  }
  static constexpr auto get_parameter (peak_frequency_tag)
  {
    return frequency_parameter (100.0, 10000.0, 880.0);
  }
  //----------------------------------------------------------------------------
  struct highshelf_frequency_tag {};
  void set (highshelf_frequency_tag, float v)
  {
    v = midi_note_to_hz (v);
    _reverb->setHighShelfFrequency ((v - 100.f) * (1.f / 9900.f));
  }
  static constexpr auto get_parameter (highshelf_frequency_tag)
  {
    return frequency_parameter (100.0, 10000.0, 2880.0);
  }
  //----------------------------------------------------------------------------
  struct lowshelf_gain_tag {};
  void set (lowshelf_gain_tag, float v)
  {
    _reverb->setLowShelfGain (rescale_eq_gain (v));
  }
  static constexpr auto get_parameter (lowshelf_gain_tag)
  {
    return float_param ("dB", 0., 18., 0.0, 0.01);
  }
  //----------------------------------------------------------------------------
  struct peak_gain_tag {};
  void set (peak_gain_tag, float v)
  {
    _reverb->setPeakGain (rescale_eq_gain (v));
  }
  static constexpr auto get_parameter (peak_gain_tag)
  {
    return float_param ("dB", 0., 18., 0.0, 0.01);
  }
  //----------------------------------------------------------------------------
  struct highshelf_gain_tag {};
  void set (highshelf_gain_tag, float v)
  {
    _reverb->setHighShelfGain (rescale_eq_gain (v));
  }
  static constexpr auto get_parameter (highshelf_gain_tag)
  {
    return float_param ("dB", 0., 18., 0.0, 0.01);
  }
  //----------------------------------------------------------------------------
  struct stereo_width_tag {};
  void set (stereo_width_tag, float v) { _reverb->setStereoWidth (v * 0.01); }
  static constexpr auto get_parameter (stereo_width_tag)
  {
    return float_param ("%", 0.0, 100.0, 40.0, 0.1);
  }
  //----------------------------------------------------------------------------
  using add_ducker::get_parameter;
  using add_ducker::set;
  using ducking_speed_tag     = add_ducker::ducking_speed_tag;
  using ducking_threshold_tag = add_ducker::ducking_threshold_tag;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    decay_time_tag,
    predelay_tag,
    predelay_sync_tag,
    lowshelf_frequency_tag,
    lowshelf_gain_tag,
    peak_frequency_tag,
    peak_gain_tag,
    highshelf_frequency_tag,
    highshelf_gain_tag,
    stereo_width_tag,
    ducking_speed_tag,
    ducking_threshold_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    _reverb.emplace (pc.get_sample_rate());
    _reverb->setStereoMode (true);
    _reverb->setWet (0.5f);
    _reverb->setDry (0.f);
    _predelay      = get_parameter (predelay_tag {}).defaultv;
    _predelay_sync = get_parameter (predelay_sync_tag {}).defaultv;

    add_ducker::reset (pc.get_sample_rate());

    mp11::mp_for_each<parameters> ([=] (auto param) {
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    add_ducker::process (
      outs,
      ins,
      samples,
      [=] (crange<T*> outs_fw, crange<T const*> ins_fw, uint samples_fw) {
        this->process_intern (outs_fw, ins_fw, samples_fw);
      });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  void process_intern (
    crange<T*>       outs,
    crange<T const*> ins,
    uint             block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    for (uint i = 0; i < block_samples; ++i) {
      float l, r;
      l = (float) ins[0][i];
      r = (float) ins[1][i];
      _reverb->process (&l, &r);
      outs[0][i] = l;
      outs[1][i] = r;
    }
  }
  //----------------------------------------------------------------------------
  void update_predelay()
  {
    // syncing in 1/16 steps
    float sync_ms = (60. / (_plugcontext->get_play_state().bpm * 16.)) * 1000.;
    float pd      = (std::floor (_predelay_sync / sync_ms) * sync_ms);
    pd += _predelay;
    pd = std::clamp (
      pd,
      get_parameter (predelay_tag {}).min, // 0
      get_parameter (predelay_tag {}).max); // 1000ms
    _reverb->setPreDelay (pd * 0.001);
  }

  //----------------------------------------------------------------------------
  static float rescale_eq_gain (float v)
  {
    constexpr float minus18dBGain = 0.12589254117941673f;
    constexpr float minus0dBGain  = 1.f;

    // reminder: linear rescaling:
    //
    //        (b - a) * (x - min)
    // f(x) = -------------------- + a
    //             (max - min)
    //
    // Where the new range is [a, b] and the old range is [min, max]

    v = db_to_gain (-v);
    v = (v - minus18dBGain) * (1.f / (minus0dBGain - minus18dBGain));
    return v;
  }

  plugin_context*                                           _plugcontext;
  float                                                     _predelay;
  float                                                     _predelay_sync;
  std::optional<::artv_dsp_pull::tal_reverb2::ReverbEngine> _reverb;
};
//------------------------------------------------------------------------------
}} // namespace artv::tal
