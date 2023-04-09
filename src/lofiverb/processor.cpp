#include <cassert>
#include <cmath>
#include <cstdio>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/dsp/own/classes/mix.hpp"
#include "artv-common/juce/effect_base.hpp"
#include "artv-common/juce/math.hpp"
#include "artv-common/juce/plugin_context.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

#include "lofiverb/parameters.hpp"

namespace artv {
// -----------------------------------------------------------------------------
// declared on editor.cpp
extern juce::AudioProcessorEditor* new_editor (
  juce::AudioProcessor&               p,
  juce::AudioProcessorValueTreeState& params,
  juce::ValueTree&                    gui_params);
// -----------------------------------------------------------------------------
class processor
  : public effect_base,
    private has_processor_params<parameters::parameters_typelist> {
public:
  //----------------------------------------------------------------------------
  processor() : effect_base {get_default_bus_properties(), make_apvts_layout()}
  {
    this->init_aptvs_references (params);
    this->setCurrentProgram (0);
  }
  //----------------------------------------------------------------------------
  juce::AudioProcessorEditor* createEditor() override
  {
    return new_editor (*this, this->params, this->gui_data);
  }
  //----------------------------------------------------------------------------
  void processBlock (juce::AudioBuffer<float>& samples, juce::MidiBuffer& midi)
    override
  {
    juce::ignoreUnused (midi);
    juce::FloatVectorOperations::disableDenormalisedNumberSupport (true);

    // read parameters rom apvts
    mp11::mp_for_each<parameters::lofiverb_parameters> ([=] (auto param) {
      using param_type = decltype (param);
      using param_tag  = typename decltype (param_type::common)::dsp_param;
      auto v           = p_refresh (param, 0);
      _fx.set (param_tag {}, v.current);
    });

    mp11::mp_for_each<parameters::mixer_parameters> ([=] (auto param) {
      using param_type = decltype (param);
      using param_tag  = typename decltype (param_type::common)::dsp_param;
      auto v           = p_refresh (param, 0);
      _mixer.set (param_tag {}, v.current);
    });

    assert (samples.getNumChannels() == 2);
    uint block_size = samples.getNumSamples();

    array2d<float, max_blocksize, 2> wetmem;
    std::array<float*, 2>            dry
      = {samples.getWritePointer (0), samples.getWritePointer (1)};
    std::array<float*, 2> wet = {wetmem[0].data(), wetmem[1].data()};

    while (block_size > 0) {
      uint                        n_spls = std::min (block_size, max_blocksize);
      std::array<float const*, 2> constdry = {dry[0], dry[1]};
      _fx.process (xspan {wet}, xspan {constdry}, n_spls);
      std::array<float const*, 4> mixerins = {dry[0], dry[1], wet[0], wet[1]};
      _mixer.process_block (xspan {dry}, xspan {mixerins}, n_spls);
      dry[0] += n_spls;
      dry[1] += n_spls;
      block_size -= n_spls;
    }
  }
  //----------------------------------------------------------------------------
  void prepareToPlay (double samplerate, int max_block_samples) override
  {
    _fx_context.reset (*this, samplerate, max_blocksize);
    _fx.reset (_fx_context);
    _mixer.reset (_fx_context);
    _mixer.skip_smoothing();
  }
  //----------------------------------------------------------------------------
#if 0
  int getNumPrograms() override
  {
    // return sizeof presets / sizeof (preset);
  }
#endif
  //----------------------------------------------------------------------------
  int getCurrentProgram() override { return _program; }
  //----------------------------------------------------------------------------
  void setCurrentProgram (int index) override
  {
    // assert (index >= 0 && index < getNumPrograms());
    _program = index;
    // TODO:
    // set_state (juce::parseXML (presets[_program].xml));
  }
  //----------------------------------------------------------------------------
  void setStateInformation (void const* src, int src_bytes) override
  {
    effect_base::setStateInformation (src, src_bytes);
  }
  //----------------------------------------------------------------------------
  const juce::String getProgramName (int index) override
  {
    assert (index >= 0 && index < getNumPrograms());
    return "TBD";
    // TODO
    // return presets[index].name;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static BusesProperties get_default_bus_properties()
  {
    BusesProperties b;
    b.addBus (true, "In", juce::AudioChannelSet::stereo(), true);
    b.addBus (false, "Out", juce::AudioChannelSet::stereo(), true);
    return b;
  }
  //----------------------------------------------------------------------------
  static constexpr uint max_blocksize = 64;
  //----------------------------------------------------------------------------
  juce_plugin_context _fx_context;
  lofiverb            _fx;
  dry_wet_mixer       _mixer;
  uint                _program {};
};

} // namespace artv
// -----------------------------------------------------------------------------

// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
  return new artv::processor {};
}
