#pragma once

#include "artv-common/misc/short_ints.hpp"

#include <juce_audio_processors/juce_audio_processors.h>
#include <utility>

namespace artv {

// processor -------------------------------------------------------------------
class effect_base : public juce::AudioProcessor {
public:
  using param_layout = juce::AudioProcessorValueTreeState::ParameterLayout;
  //----------------------------------------------------------------------------
  effect_base (
    juce::AudioProcessor::BusesProperties const& bp,
    param_layout&&                               pl)
    : AudioProcessor (bp)
    , params {*this, nullptr, juce::Identifier ("params"), std::forward<param_layout> (pl)}
    , gui_data {"gui_data"}
  {}
  //----------------------------------------------------------------------------
  effect_base (effect_base const&) = delete;
  effect_base& operator= (effect_base const&) = delete;

  ~effect_base() override {}
  //----------------------------------------------------------------------------
  void prepareToPlay (double sampleRate, int samplesPerBlock) override
  {
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
    juce::ignoreUnused (sampleRate, samplesPerBlock);
  }
  //----------------------------------------------------------------------------
  void releaseResources() override
  {
    // When playback stops, you can use this as an opportunity to free up
    // any spare memory, etc.
  }
  //----------------------------------------------------------------------------
  bool isBusesLayoutSupported (const BusesLayout& layouts) const override
  {
#if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
#else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    auto mono   = juce::AudioChannelSet::mono();
    auto stereo = juce::AudioChannelSet::stereo();
    auto outset = layouts.getMainOutputChannelSet();
    if (outset != mono && outset != stereo) {
      return false;
    }
// This checks if the input layout matches the output layout
#if !JucePlugin_IsSynth
    if (outset != layouts.getMainInputChannelSet()) {
      return false;
    }
#endif
    return true;
#endif
  }
  //----------------------------------------------------------------------------
  void processBlock (juce::AudioBuffer<float>& samples, juce::MidiBuffer& midi)
    override
  {
    juce::ignoreUnused (samples, midi);
#if 0 // remove this


        juce::ScopedNoDenormals noDenormals;
        auto in_count  = getTotalNumInputChannels();
        auto out_count = getTotalNumOutputChannels();

        // In case we have more outputs than inputs, this code clears any output
        // channels that didn't contain input data, (because these aren't
        // guaranteed to be empty - they may contain garbage).
        // This is here to avoid people getting screaming feedback
        // when they first compile a plugin, but obviously you don't need to
        // keep this code if your algorithm always overwrites all the output
        // channels.
        for (uint i = in_count; i < out_count; ++i)
            samples.clear (i, 0, samples.getNumSamples());

        // This is the place where you'd normally do the guts of your plugin's
        // audio processing...
        // Make sure to reset the state if your inner loop is processing
        // the samples and the outer loop is handling the channels.
        // Alternatively, you can process the samples with the channels
        // interleaved by keeping the same state.
        for (int channel = 0; channel < in_count; ++channel)
        {
            auto* channelData = samples.getWritePointer (channel);
            juce::ignoreUnused (channelData);
            // ..do something to the data...
        }
#endif
  }
  //----------------------------------------------------------------------------
  bool hasEditor() const override { return true; }
  //----------------------------------------------------------------------------
  const juce::String getName() const override { return JucePlugin_Name; }
  //----------------------------------------------------------------------------
  bool acceptsMidi() const override
  {
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
  }
  //----------------------------------------------------------------------------
  bool producesMidi() const override
  {
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
  }
  //----------------------------------------------------------------------------
  bool isMidiEffect() const override
  {
#if JucePlugin_IsMidiEffect
    return true;
#else
    return false;
#endif
  }
  //----------------------------------------------------------------------------
  double getTailLengthSeconds() const override { return 0.0; }
  //----------------------------------------------------------------------------
  int getNumPrograms() override
  {
    // NB: some hosts don't cope very well if you tell them there are 0
    // programs, so this should be at least 1, even if you're not really
    // implementing programs.
    return 1;
  }
  //----------------------------------------------------------------------------
  int getCurrentProgram() override { return 0; }
  //----------------------------------------------------------------------------
  void setCurrentProgram (int index) override { juce::ignoreUnused (index); }
  //----------------------------------------------------------------------------
  const juce::String getProgramName (int index) override
  {
    juce::ignoreUnused (index);
    return {};
  }
  //----------------------------------------------------------------------------
  void changeProgramName (int index, const juce::String& newName) override
  {
    juce::ignoreUnused (index, newName);
  }
  //----------------------------------------------------------------------------
  void getStateInformation (juce::MemoryBlock& dst) override
  {
    using xml_ptr = std::unique_ptr<juce::XmlElement>;
    xml_ptr xml {params.copyState().createXml()};
    xml_ptr gui_xml {gui_data.createXml()};
    if (gui_xml) {
      auto edit = xml->createNewChildElement (gui_data.getType());
      *edit     = *gui_xml;
    }
#if !defined(NDEBUG) && 0
    printf ("Save preset\n%s\n", (char const*) xml->toString().toUTF8());
#endif
    copyXmlToBinary (*xml, dst);
  }
  //----------------------------------------------------------------------------
  void setStateInformation (void const* src, int src_bytes) override
  {
    set_state (getXmlFromBinary (src, src_bytes));
  }
  //----------------------------------------------------------------------------
  void set_state (std::unique_ptr<juce::XmlElement> xml)
  {
    if (!xml) {
      return;
    }
#if !defined(NDEBUG) && 0
    printf ("Set preset\n%s\n", (char const*) xml->toString().toUTF8());
#endif
    juce::XmlElement* gui = xml->getChildByName (gui_data.getType());
    if (gui) {
      // update expected values
      auto gui_new = juce::ValueTree::fromXml (*gui);
      for (uint i = 0; i < gui_data.getNumProperties(); ++i) {
        auto key  = gui_data.getPropertyName (i);
        auto prop = gui_new.getPropertyPointer (key);
        if (prop) {
          gui_data.setProperty (key, *prop, nullptr);
          gui_new.removeProperty (key, nullptr);
        }
        else {
          // not present means empty/default.
          gui_data.setProperty (key, "", nullptr);
        }
      }
      // keep unknown values on the preset, so an old version doesn't destroy
      // new presets.
      for (uint i = 0; i < gui_new.getNumProperties(); ++i) {
        auto key = gui_new.getPropertyName (i);
        gui_data.setProperty (key, gui_new.getProperty (key), nullptr);
      }
      xml->removeChildElement (gui, true);
    }
    else {
      // clear everything
      for (uint i = 0; i < gui_data.getNumProperties(); ++i) {
        auto name = gui_data.getPropertyName (i);
        gui_data.setProperty (name, "", nullptr);
      }
    }
    if (xml->hasTagName (params.state.getType())) {
      params.replaceState (juce::ValueTree::fromXml (*xml));
    }
    preset_was_loaded (std::move (xml));
  }
  //----------------------------------------------------------------------------
  virtual void preset_was_loaded (std::unique_ptr<juce::XmlElement>) {}
  //----------------------------------------------------------------------------
  juce::AudioProcessorValueTreeState params;
  juce::ValueTree                    gui_data;
};

} // namespace artv
