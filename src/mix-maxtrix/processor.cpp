#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdint.h>
#include <utility>
#include <variant>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include <ff_meters.h>

#include "artv-common/dsp/own/classes/mix.hpp"
#include "artv-common/juce/effect_base.hpp"
#include "artv-common/juce/math.hpp"
#include "artv-common/juce/parameters.hpp"
#include "artv-common/juce/plugin_context.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/delay_compensation_buffers.hpp"
#include "artv-common/misc/short_ints.hpp"

#include "mix-maxtrix/buffers.hpp"
#include "mix-maxtrix/order_and_buffering.hpp"
#include "mix-maxtrix/parameters.hpp"
#include "mix-maxtrix/presets.hpp"

namespace artv {
// -----------------------------------------------------------------------------
class mix_maxtrix_fx_context : public juce_plugin_context {
public:
  void reset (juce::AudioProcessor& proc, uint samplerate, uint max_block_size)
  {
    juce_plugin_context::reset (proc, samplerate, max_block_size);
    memset (&fx_latencies, 0, sizeof fx_latencies);
    fx_latencies_changed = false;
    _processing_chnl_idx = 0;
  }
  // ---------------------------------------------------------------------------
  void set_delay_compensation (uint samples) override
  {
    fx_latencies_changed |= fx_latencies[_processing_chnl_idx] != samples;
    fx_latencies[_processing_chnl_idx] = samples;
  }
  // ---------------------------------------------------------------------------
  void set_processing_channel (uint chnl) { _processing_chnl_idx = chnl; }
  // ---------------------------------------------------------------------------
  void on_fx_disable (uint chnl)
  {
    fx_latencies_changed |= fx_latencies[chnl] != 0;
    fx_latencies[chnl] = 0;
  }
  // ---------------------------------------------------------------------------
  order_and_buffering<parameters::n_stereo_busses>::bus_latency_arr
       fx_latencies;
  bool fx_latencies_changed;
  // ---------------------------------------------------------------------------
private:
  uint _processing_chnl_idx = 0;
};
// -----------------------------------------------------------------------------
template <class T>
using get_dsp_type = stereo_processor_adapt<
  typename decltype (mp11::mp_front<T>::common)::dsp_type>;

using dsp_class_typelist
  = mp11::mp_transform<get_dsp_type, parameters::all_fx_typelists>;

using dsp_variant_class_typelist
  = mp11::mp_push_front<dsp_class_typelist, std::monostate>;

using dsp_variant = mp11::mp_rename<dsp_variant_class_typelist, std::variant>;

constexpr size_t dsp_variant_size = mp11::mp_size<dsp_variant>::value;

using indexed_all_fx_params_typelists
  = to_typelist_with_indexes<parameters::all_fx_typelists>;
// -----------------------------------------------------------------------------
// declared on editor.cpp
juce::AudioProcessorEditor* new_editor (
  juce::AudioProcessor&               p,
  juce::AudioProcessorValueTreeState& params,
  crange<foleys::LevelMeterSource>    meter_srcs);
// -----------------------------------------------------------------------------
class processor
  : public effect_base,
    private has_processor_params<parameters::parameters_typelist> {
public:
  //----------------------------------------------------------------------------
  processor()
    : effect_base {get_default_bus_properties(), make_apvts_layout()}
    , _fsamples {
        parameters::n_stereo_busses * 2,
        (parameters::crossover_n_bands - 1) * 2}
  {
    this->init_aptvs_references (params);
    for (auto& meter_src : _meters) {
      meter_src.setSuspended (true);
    }
    // adding version
    params.state.setProperty (
      juce::Identifier ("plugin_version"), VERSION_INT, nullptr);
    this->setCurrentProgram (0);

    // register listeners for parameter changes.
    mp11::mp_for_each<parameters::routing_controls_typelist> ([=] (auto param) {
      for (char const* id : param.juce_ids) {
        params.addParameterListener (id, &_routing_refreshed);
      }
    });
    mp11::mp_for_each<parameters::non_routing_controls_typelist> (
      [=] (auto param) {
        for (char const* id : param.juce_ids) {
          params.addParameterListener (id, &_non_routing_refreshed);
        }
      });
    mp11::mp_for_each<parameters::all_nonfx_sliders_typelist> (
      [=] (auto param) {
        for (char const* id : param.juce_ids) {
          params.addParameterListener (id, &_non_fx_slider_refreshed);
        }
      });
  }
  //----------------------------------------------------------------------------
  juce::AudioProcessorEditor* createEditor() override
  {
    return new_editor (*this, this->params, this->_meters);
  }
  //----------------------------------------------------------------------------
  void processBlock (juce::AudioBuffer<float>& samples, juce::MidiBuffer& midi)
    override
  {
    juce::ignoreUnused (midi);
    process_block (samples);
  }
  //----------------------------------------------------------------------------
  void prepareToPlay (double samplerate, int max_block_samples) override
  {
    _fx_context.reset (*this, samplerate, fx_blocksize);
    _crossv_type = crossv_off;

    for (uint i = 0; i < parameters::n_stereo_busses; ++i) {
      _mixer[i].reset (_fx_context);
      _mixer[i].skip_smoothing();
      _meters[i].resize (2, samplerate * 0.1 / max_block_samples);

      mp11::mp_with_index<dsp_variant_size> (
        _fx_dsp[i].index(), [=] (auto idx) {
          auto& fx       = std::get<idx> (_fx_dsp[i]);
          using dsp_type = std::remove_reference_t<decltype (fx)>;
          if constexpr (!std::is_same<dsp_type, std::monostate>::value) {
            fx.reset (_fx_context);
          }
        });
    }
    // for "io_recompute" to forcefully run on first block we set an invalid
    // value on one variable affecting the routing, so when it is refreshed it
    // // will always have a change. This hack is done to avoid code bloating
    // this function.
    p_get (parameters::in_selection {})[0] = 1 << parameters::n_stereo_busses;
    // force reloading of all parameters
    _routing_refreshed.parameterChanged ("", 0.f);
    _non_routing_refreshed.parameterChanged ("", 0.f);
    _non_fx_slider_refreshed.parameterChanged ("", 0.f);
  }
  //----------------------------------------------------------------------------
  bool isBusesLayoutSupported (const BusesLayout& buses) const override
  {
    using namespace juce;

    if (buses.inputBuses.size() > parameters::n_stereo_busses) {
      return false;
    }
    if (buses.outputBuses.size() > parameters::n_stereo_busses) {
      return false;
    }

    AudioChannelSet stereo   = AudioChannelSet::stereo();
    AudioChannelSet disabled = AudioChannelSet::disabled();

    for (int bus = 0; bus < parameters::n_stereo_busses; ++bus) {
      // only stereo or disabled buses.
      auto const& inbus  = buses.getChannelSet (true, bus);
      auto const& outbus = buses.getChannelSet (false, bus);

      if (inbus != stereo || inbus != disabled) {
        return false;
      }
      if (outbus != stereo || outbus != disabled) {
        return false;
      }
    }
    return true;
  }
  //----------------------------------------------------------------------------
  int getNumPrograms() override { return sizeof presets / sizeof (preset); }
  //----------------------------------------------------------------------------
  int getCurrentProgram() override { return _program; }
  //----------------------------------------------------------------------------
  void setCurrentProgram (int index) override
  {
    assert (index >= 0 && index < getNumPrograms());
    _program = index;
    set_apvts_state (juce::parseXML (presets[_program].xml));
  }
  //----------------------------------------------------------------------------
  const juce::String getProgramName (int index) override
  {
    assert (index >= 0 && index < getNumPrograms());
    return presets[index].name;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  buffers<float>& get_buffers (float) { return _fsamples; }
  //----------------------------------------------------------------------------
  auto& get_mix_buffer (float) { return _fx_mix_buffer; }
  //----------------------------------------------------------------------------
  static BusesProperties get_default_bus_properties()
  {
    BusesProperties b;
    for (uint i = 1; i <= parameters::n_stereo_busses; ++i) {
      auto n = juce::String (i);
      b.addBus (true, "In-" + n, juce::AudioChannelSet::stereo(), true);
      b.addBus (false, "Out-" + n, juce::AudioChannelSet::stereo(), true);
    }
    return b;
  }
  //----------------------------------------------------------------------------
  struct mutesolo_state {
    bool changed;
    u64  bits;
  };
  //----------------------------------------------------------------------------
  void io_recompute (
    u64                                                  mute_solo_bits,
    std::array<uint, parameters::crossover_n_bands - 1>& crossv_outs)
  {
    auto& ins_arr      = p_get (parameters::in_selection {});
    auto& outs_arr     = p_get (parameters::out_selection {});
    auto& mixsends_arr = p_get (parameters::mixer_sends {});
    auto& routing_arr  = p_get (parameters::routing {});

    u64 ins  = 0;
    u64 outs = 0;

    assert (ins_arr.size() == outs_arr.size());
    for (uint i = 0; i < ins_arr.size(); ++i) {
      ins |= (((u64) ins_arr[i]) << (i * parameters::n_stereo_busses));
      outs |= (((u64) outs_arr[i]) << (i * parameters::n_stereo_busses));
    }

    // if the channel 0, which is the only one that can have  the crossover, has
    // no inputs, then the crossover is off.
    _crossv_type     = (ins_arr[0] == 0) ? crossv_off : _crossv_type;
    uint inhibit_opt = 0;
    if (_crossv_type != crossv_off) {
      // This was never thought out to be a multifx, so it was never was
      // designed to have multiple ins (sidechain) and outs per fx. It is
      // perfectly doable backend-wise but not on the GUI. The crossover is
      // added by implementing 3 send channels by inhibiting reordering
      // optimizations on the affected channels and manually setting the 3
      // additional crossover's buffers as inputs to the channel. See the
      // conditional depending on "_crossv_type" below.
      for (uint i = 0; i < crossv_outs.size(); ++i) {
        inhibit_opt |= 1 << crossv_outs[i];
      }
    }

    uint routing = (uint) routing_arr[0];

    auto latency_prev = _io.plugin_latency;
    _io.recompute (
      ins,
      outs,
      mute_solo_bits,
      mixsends_arr[0],
      inhibit_opt,
      routing,
      _fx_context.fx_latencies);

    if (_crossv_type != crossv_off) {
      // buffers -9,-10 and -11 are reserved for crossover outputs.
      uint limit = parameters::n_parallel_buses (routing);
      for (uint i = 0; i < crossv_outs.size(); ++i) {
        if (crossv_outs[i] < limit) {
          // only sending to the first group of parallel busses.
          _io.in[crossv_outs[i]][_io.n_busses + i]
            = -((int) _io.n_busses) - 1 - i;
        }
      }
    }

    this->setLatencySamples (_io.plugin_latency);
    _fx_context.fx_latencies_changed = false;
    reset_latency_buffers();
    clear_crossv_buffers();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block (juce::AudioBuffer<T>& passed_buffers)
  {
    using value_type = T;
    // I think that this is coarse grained/greedy, but
    // "juce::ScopedNoDenormals" doesn't work (Reaper).
    juce::FloatVectorOperations::disableDenormalisedNumberSupport (true);

    auto& buffers = get_buffers (T {});
    buffers.on_block (passed_buffers);

    bool                                          routing_change = false;
    mutesolo_state                                mutesolo {false, 0};
    std::array<bool, parameters::n_stereo_busses> fx_type_changed = {};

    // move (refresh) all non-fx parameters from JUCE's atomic_ptrs to local
    // member variables, so we read the same value for the whole duration on
    // this block. There is a substantial amount of parameters. These are
    // refreshed in declaration order to avoid weird memory access orders.
    if (unlikely (_routing_refreshed.get_and_clear())) {
      p_refresh_many (
        parameters::routing_controls_typelist {},
        [&] (auto key, uint i, auto val) { routing_change |= val.changed(); });
    }

    bool crossv_change = false;
    if (unlikely (
          _non_routing_refreshed.get_and_clear() || routing_change
          || _fx_context.fx_latencies_changed)) {
      // routing or latency chanes have to re-read the mute solo state because
      // io recompute will need the mute_solo bits.
      p_refresh_many (
        parameters::non_routing_controls_typelist {},
        [&] (auto key, uint i, auto val) {
          if constexpr (std::is_same_v<decltype (key), parameters::fx_type>) {
            constexpr uint crossv_id
              = mp11::mp_find<
                  parameters::all_fx_typelists,
                  parameters::lr_crossv_params>::value
              + 1;
            constexpr uint wonky_crossv_id
              = mp11::mp_find<
                  parameters::all_fx_typelists,
                  parameters::wonky_crossv_params>::value
              + 1;

            fx_type_changed[i] = val.changed();
            if (i != 0) {
              return; // crossover only on channel 0
            }
            if (val.current == crossv_id) {
              _crossv_type = crossv_normal;
            }
            else if (val.current == wonky_crossv_id) {
              _crossv_type = crossv_wonky;
            }
            crossv_change = fx_type_changed[i]
              && (_crossv_type != crossv_off || val.previous == crossv_id
                  || val.previous == wonky_crossv_id);
          }
          if constexpr (std::is_same_v<decltype (key), parameters::mute_solo>) {
            mutesolo.changed |= val.changed();
            mutesolo.bits |= ((u64) val.current)
              << (i * 2); // 2 mute solo bits.
          }
        });
    }

    std::array<uint, parameters::crossover_n_bands - 1> crossv_outs {};
    if (unlikely (_crossv_type == crossv_normal)) {
      using band_params = mp11::mp_list<
        parameters::lr_crossv_band1_out,
        parameters::lr_crossv_band2_out,
        parameters::lr_crossv_band3_out>;

      mp_foreach_idx (band_params {}, [&] (auto index, auto param) {
        auto value = p_refresh (param, 0);
        // notice: +1 because the crossover can't output on channel 0
        crossv_outs[index] = value.current + 1;
        crossv_change |= value.changed();
      });
    }
    else if (unlikely (_crossv_type == crossv_wonky)) {
      using band_params = mp11::mp_list<
        parameters::wonky_crossv_band1_out,
        parameters::wonky_crossv_band2_out,
        parameters::wonky_crossv_band3_out>;

      mp_foreach_idx (band_params {}, [&] (auto index, auto param) {
        auto value = p_refresh (param, 0);
        // notice: +1 because the crossover can't output on channel 0
        crossv_outs[index] = value.current + 1;
        crossv_change |= value.changed();
      });
    }

    if (unlikely (_non_fx_slider_refreshed.get_and_clear())) {
      p_refresh_many (parameters::all_nonfx_sliders_typelist {});
    }

    if (unlikely (
          routing_change || _fx_context.fx_latencies_changed
          || crossv_change)) {
      // not delaying mute and solo...
      io_recompute (mutesolo.bits, crossv_outs);
      mutesolo.changed = false;
    }

    channels_mute_solo<parameters::n_stereo_busses> mixer_crossfade;

    if (likely (!mutesolo.changed)) {
      mixer_crossfade = _io.mute_solo;
    }
    else {
      mixer_crossfade.reset (mutesolo.bits);
    }

    // swap channels
    for (uint i = 0; i < _io.swaps.size(); ++i) {
      if (_io.swaps[i] != 0) {
        buffers.swap (i + 1, _io.swaps[i]);
      }
    }

    std::array<bool, parameters::n_stereo_busses> zeroed {};

    uint rmode     = p_get (parameters::routing {})[0];
    uint n_gbusses = parameters::n_parallel_buses (rmode);

    for (uint bus_beg = 0; bus_beg < _io.order.size(); bus_beg += n_gbusses) {

      auto bus_end = bus_beg + n_gbusses;
      // process mix
      for (uint i = bus_beg; i < bus_end; ++i) {
        uint chnl = _io.order[i];

        if (unlikely (!refresh_fx_params_and_mix_inputs (
              chnl, buffers, fx_type_changed[chnl]))) {
          zeroed[chnl] = true;
          continue; // content is 0.s (unconnected, etc.)
        }
        process_mix<T> (chnl, !mixer_crossfade.is_channel_active (chnl));
        process_modifiers<T> (chnl);
      }

      for (uint i = bus_beg; i < bus_end; ++i) {
        // delay compensation
        if (_io.pre_output_mix_delay_samples[i] != 0) {
          dly_compensate_pre_output_mix (i);
        }
        // metering
        uint buf_idx = buffers.get_first_audiobuffer_chnl (_io.mix[i], 2);
        _meters[i].measureBlock (
          buffers.get_audiobuffer (_io.mix[i]), buf_idx, 0, 2);
      }

      // output mix
      for (uint chnl : _io.order) {
        typename decltype (_io.out)::value_type ins {};

        for (uint i = bus_beg; i < bus_end; ++i) {
          ins[i] = zeroed[i] ? 0 : _io.out[chnl][i];
        }
        buffers.template mix<2, s8> (
          chnl + 1, make_crange (ins)); // chnls are 1 indexed
      }
    }

    if (unlikely (mutesolo.changed)) {
      // delay IO recomputation until we have crossfaded a mute. This will
      // cause a 1 audio buffer delay when unmuting + 1 audio buffer crossfade
      // when unmuting, but it is the only sensible way to do when having a
      // dynamic buffer ordering optimization.
      io_recompute (mutesolo.bits, crossv_outs);
    }
  }
  //----------------------------------------------------------------------------
  // these two seemingly unrelated things are done in place to avoid
  // generating unnecessary compile time switches on the fast-path.
  template <class T>
  bool refresh_fx_params_and_mix_inputs (
    uint        channel,
    buffers<T>& buffers,
    bool        fx_type_changed)
  {
    // + 1 for the mixer 2 mixer sends
    std::array<stereo_summing_processor*, parameters::console_n_elems>
      processors {};

    auto           fx_type = p_get (parameters::fx_type {})[channel];
    constexpr auto fx_val  = parameters::fx_type {};

    if (unlikely (fx_type_changed)) {
      if (fx_type == 0 || fx_type >= fx_val.type.choices.size()) {
        _fx_dsp[channel] = std::monostate {};
        _fx_context.on_fx_disable (channel);
      }
    }

    // TODO: work might need to be done here to ensure no code generation
    // bloat.
    // clang-format off
    mp11::mp_for_each<indexed_all_fx_params_typelists> (
        [=, &processors] (auto param_tlist_widx) {

      using param_tlist_widx_t = decltype (param_tlist_widx);
      auto fx_idx              = type_with_index_index<param_tlist_widx_t> {};
      using param_tlist        = type_with_index_type<param_tlist_widx_t>;
      using dsp_class          = get_dsp_type<param_tlist>;
      using stereo_processor_type = typename dsp_class::stereo_processor_type;

      if ((fx_idx.value + 1) != fx_type) {
        // +1 because the first entry is "none"/std::monostate
        return; // not this fx type. Keep the linear search
      }

      // setting the processing channel, so latencies get updated.
      _fx_context.set_processing_channel (channel);

      static constexpr bool uses_oversampled =
        is_same_template_v<oversampled, stereo_processor_type>;

      // refresh FX type.
      if (unlikely (fx_type_changed)) {
        _fx_dsp[channel] = dsp_class{};
        auto& fx = std::get<dsp_class> (_fx_dsp[channel]);
        fx.reset (_fx_context);
      }

      // refresh FX parameters.
      auto& fx = std::get<dsp_class> (_fx_dsp[channel]);

      if constexpr (uses_oversampled) {
        // run the oversampling parameter first (if any)
        //
        // The "oversampled" class fully resets the FX on samplerate changes.
        // This way all the parameters will be passed afterwards. The next
        // loop of parameter setting sets the oversampling again, but
        // "oversampled" just acts on changes, so it is OK to not handle the
        // case.
        mp11::mp_for_each<param_tlist> ([=, &fx] (auto param) {
          using paramtype = decltype (param);
          using dsp_param = typename decltype (paramtype::common)::dsp_param;
          static constexpr bool is_oversampling_param =
            std::is_same_v<dsp_param, oversampled_amount_tag>;

          if constexpr (is_oversampling_param) {
            auto v = p_refresh (param, channel);
            fx.set (dsp_param {}, v.current);
          }
        });
      }

      mp11::mp_for_each<param_tlist> ([=, &fx] (auto param) {
        using paramtype = decltype (param);
        using dsp_param = typename decltype (paramtype::common)::dsp_param;

        auto v = p_refresh (param, channel);
        fx.set (dsp_param {}, v.current);
      });

      // Get input summing instances on console dsp types.
      if constexpr (dsp_class::dsp_type == dsp_types::console) {
        for (uint i = 0; i < processors.size(); ++i) {
          processors[i] = &fx.get_channel (i);
        }
      }
    });
    // clang-format on
    bool ret;
    if (processors[0] == nullptr) {
      ret = buffers.template mix<2, s8> (
        _io.mix[channel], make_crange (_io.in[channel]));

      if (_io.post_input_mix_delay_samples[channel] != 0) {
        dly_compensate_post_input_mix (channel);
      }
      if (_io.receives[channel] != 0) {
        ret |= buffers.template mix<2, s8> (
          _io.mix[channel], make_crange (_io.receives[channel]), ret);
      }
    }
    else {
      ret = buffers.template mix<2, s8> (
        _io.mix[channel],
        make_crange (_io.in[channel]),
        make_crange (processors));

      if (_io.post_input_mix_delay_samples[channel] != 0) {
        dly_compensate_post_input_mix (channel);
      }
      if (_io.receives[channel] != 0) {
        ret |= buffers.template mix<2, s8> (
          _io.mix[channel],
          make_crange (_io.receives[channel]),
          make_crange (processors),
          ret);
      }
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  void mixer_parameters_update (uint channel, bool will_mute)
  {
    float wet_pan = p_get (parameters::wet_pan {})[channel];
    float dry_pan = p_get (parameters::dry_pan {})[channel];

    float wet_ms_bal = p_get (parameters::wet_balance {})[channel];
    float dry_ms_bal = p_get (parameters::dry_balance {})[channel];

    float pan    = p_get (parameters::pan {})[channel];
    float vol    = p_get (parameters::volume {})[channel];
    float gvol   = p_get (parameters::global_volume {})[0];
    int   modifs = p_get (parameters::channel_modifs {})[channel];
    float mix    = p_get (parameters::fx_mix {})[channel];

    vol  = db_to_gain (vol);
    gvol = db_to_gain (gvol);

    using dwm = dry_wet_mixer;
    _mixer[channel].set (dwm::gain_tag {}, !will_mute ? (vol * gvol) : 0.f);
    _mixer[channel].set (dwm::dry_wet_ratio_tag {}, mix);
    _mixer[channel].set (dwm::pan_tag {}, pan);
    _mixer[channel].set (dwm::dry_pan_tag {}, dry_pan);
    _mixer[channel].set (dwm::wet_pan_tag {}, wet_pan);
    _mixer[channel].set (dwm::dry_ms_ratio_tag {}, dry_ms_bal);
    _mixer[channel].set (dwm::wet_ms_ratio_tag {}, wet_ms_bal);
    _mixer[channel].set (dwm::phase_inv_l_tag {}, modifs & 1);
    _mixer[channel].set (dwm::phase_inv_r_tag {}, modifs & 1);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_mix (uint channel, bool will_mute)
  {
    mixer_parameters_update (channel, will_mute);

    // processing buffers
    auto&          buffers = get_buffers (T {});
    auto           mix_bus = buffers.get_write_ptrs (_io.mix[channel]);
    auto           samples = buffers.sample_count();
    auto           fx_type = p_get (parameters::fx_type {})[channel];
    constexpr auto fx_val  = parameters::fx_type {};

    if (fx_type == 0 || fx_type >= fx_val.type.choices.size()) {
      // Fast path: Skipping fx and dry/wet mixing.
      auto c_mix_bus = array_const_cast<T const*> (mix_bus);
      _mixer[channel].process_block<T> (mix_bus, c_mix_bus, samples);
      return;
    }

    // FX enabled. Fx processor lookup.
    stereo_processor<T>* fx = nullptr;
    mp11::mp_for_each<indexed_all_fx_params_typelists> (
      [=, &fx] (auto param_tlist_widx) {
        using param_tlist_widx_t = decltype (param_tlist_widx);
        auto fx_idx              = type_with_index_index<param_tlist_widx_t> {};
        using param_tlist        = type_with_index_type<param_tlist_widx_t>;
        using dsp_class          = get_dsp_type<param_tlist>;

        if ((fx_idx.value + 1) != fx_type) {
          // +1 because the first entry is "none"/std::monostate
          return;
        }
        fx = &std::get<dsp_class> (_fx_dsp[channel]);
      });

    // FX process and mix.
    assert (fx);
    auto& mix_buffer = get_mix_buffer (T {});

    constexpr auto crossv_head = parameters::n_stereo_busses * 2;
    // getting the buffer for the negative indexes, where we store the mix
    // buffers not provided by Juce and the Crossover buffers
    auto& internal_audiobuff = buffers.get_audiobuffer (-1);
    static_assert (parameters::crossover_n_bands == 4, "Fix this");
    auto out_ptrs = make_array<T*> (
      mix_buffer[0].data(),
      mix_buffer[1].data(),
      internal_audiobuff.getWritePointer (crossv_head),
      internal_audiobuff.getWritePointer (crossv_head + 1),
      internal_audiobuff.getWritePointer (crossv_head + 2),
      internal_audiobuff.getWritePointer (crossv_head + 3),
      internal_audiobuff.getWritePointer (crossv_head + 4),
      internal_audiobuff.getWritePointer (crossv_head + 5));

    for (uint i = 0; i < samples;) {
      uint bsz = std::min<uint> (samples - i, fx_blocksize);

      auto c_mix_bus = array_const_cast<T const*> (mix_bus);

      fx->process (crange<T*> {out_ptrs}, crange<T const*> {c_mix_bus}, bsz);

      if (_io.fx_delay_samples[channel] != 0) {
        dly_compensate_fx_dry (channel, mix_bus, bsz);
      }

      auto mixer_in = make_array<T const*> (
        mix_bus[0], mix_bus[1], out_ptrs[0], out_ptrs[1]);

      _mixer[channel].process_block<T> (mix_bus, mixer_in, bsz);

      i += bsz;
      // advance non-blockwise buffers.
      mix_bus[0] += bsz;
      mix_bus[1] += bsz;
      for (uint j = 2; j < out_ptrs.size(); ++j) {
        out_ptrs[j] += bsz;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_modifiers (uint channel)
  {
    auto& buffers = get_buffers (T {});
    auto  buses   = buffers.to_mono_buses (_io.mix[channel]);
    // to avoid bloat
    std::array<s8, 2> buses_s8 = {(s8) buses[0], (s8) buses[1]};

    auto modifs = (uint) p_get (parameters::channel_modifs {})[channel];
    // Phase handled on the volumes by negating the gain, bit shift +1
    modifs >>= 1;

    // Ph, L, R, RL
    switch (modifs) {
    case 0:
      return;
    case 1:
      // copy L to R
      buffers.template mix<1, s8> (buses_s8[1], make_crange (&buses_s8[0], 1));
      return;
    case 2:
      // copy R to L
      buffers.template mix<1, s8> (buses_s8[0], make_crange (&buses_s8[1], 1));
      return;
    case 4:
      buffers.template swap<1> (buses_s8[0], buses_s8[1]);
      return;
    default:
      assert (false);
      return;
    }
  }
  //----------------------------------------------------------------------------
  void reset_latency_buffers()
  {
    // placing on the memory buffer in the most likely processing order:
    //
    // for all channels (post input mix + fx)
    constexpr auto nbuses = parameters::n_stereo_busses;

    std::array<uint, nbuses * 3> lat;

    for (uint i = 0; i < nbuses; ++i) {
      uint off              = i * 2;
      lat[off]              = _io.post_input_mix_delay_samples[i];
      lat[off + 1]          = _io.fx_delay_samples[i];
      lat[(nbuses * 2) + i] = _io.pre_output_mix_delay_samples[i];
    }

    // creates "lat.size() * 2" buffers for L and R.
    _dly_comp_buffers.reset<uint> (lat);
  }
  //----------------------------------------------------------------------------
  void dly_compensate_post_input_mix (uint mix_chnl)
  {
    uint samples  = _fsamples.sample_count();
    auto chnlptrs = _fsamples.get_write_ptrs (_io.mix[mix_chnl]);
    uint id_l     = mix_chnl * (2 + 2);
    uint id_r     = id_l + 1;
    _dly_comp_buffers.compensate (id_l, make_crange (chnlptrs[0], samples));
    _dly_comp_buffers.compensate (id_r, make_crange (chnlptrs[1], samples));
  }
  //----------------------------------------------------------------------------
  void dly_compensate_fx_dry (
    uint                  mix_chnl,
    std::array<float*, 2> buff,
    uint                  n_samples)
  {
    uint id_l = (mix_chnl * (2 + 2)) + 2;
    uint id_r = id_l + 1;
    _dly_comp_buffers.compensate (id_l, make_crange (buff[0], n_samples));
    _dly_comp_buffers.compensate (id_r, make_crange (buff[1], n_samples));
  }
  //----------------------------------------------------------------------------
  void dly_compensate_pre_output_mix (uint mix_chnl)
  {
    uint samples  = _fsamples.sample_count();
    auto chnlptrs = _fsamples.get_write_ptrs (_io.mix[mix_chnl]);
    uint id_l     = (parameters::n_stereo_busses * (2 + 2)) + (mix_chnl * 2);
    uint id_r     = id_l + 1;
    _dly_comp_buffers.compensate (id_l, make_crange (chnlptrs[0], samples));
    _dly_comp_buffers.compensate (id_r, make_crange (chnlptrs[1], samples));
  }
  //----------------------------------------------------------------------------
  void clear_crossv_buffers()
  {
    auto& b = get_buffers (float {});
    for (uint i = 0; i < _io.n_external_ins; ++i) {
      b.clear (-(_io.n_busses + 1 + i));
    }
  }
  //----------------------------------------------------------------------------
  struct param_change_counter
    : public juce::AudioProcessorValueTreeState::Listener {
    virtual ~param_change_counter() {};
    virtual void parameterChanged (juce::String const&, float) override
    {
      _count.fetch_add (1, std::memory_order_release);
    }

    uint get_and_clear()
    {
      return _count.exchange (0, std::memory_order_acquire);
    }

  private:
    std::atomic<uint> _count {0};
  };
  //----------------------------------------------------------------------------
  using io_engine = order_and_buffering<
    parameters::n_stereo_busses,
    parameters::crossover_n_bands - 1>;

  mix_maxtrix_fx_context                                 _fx_context;
  buffers<float>                                         _fsamples;
  io_engine                                              _io;
  std::array<dry_wet_mixer, parameters::n_stereo_busses> _mixer;

  std::array<dsp_variant, parameters::n_stereo_busses> _fx_dsp;

  static constexpr uint fx_blocksize = 128;
  alignas (32) std::array<std::array<float, fx_blocksize>, 2> _fx_mix_buffer;
  delay_compensation_buffers<float, 2> _dly_comp_buffers;
  uint                                 _program = 0;

  enum crossv_type {
    crossv_off,
    crossv_normal,
    crossv_wonky,
  };
  uint _crossv_type;

  alignas (128)
    std::array<foleys::LevelMeterSource, parameters::n_stereo_busses> _meters;

  alignas (128) param_change_counter _routing_refreshed;
  param_change_counter _non_routing_refreshed;
  param_change_counter _non_fx_slider_refreshed;
};

} // namespace artv
// -----------------------------------------------------------------------------

// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  return new artv::processor();
}
