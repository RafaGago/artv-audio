#include <optional>

#include <cerrno>
#include <fstream>
#include <streambuf>

#include <gtest/gtest.h>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/misc/short_ints.hpp"

#if 0
#include <cstdio>
#define dbg(...) printf (__VA_ARGS__)
#else
#define dbg(...)
#endif

extern juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter();

namespace artv {

std::string get_file_contents (char const* filename)
{
  std::ifstream in (filename, std::ios::in | std::ios::binary);
  if (!in) {
    throw (errno);
  }
  return (std::string (
    (std::istreambuf_iterator<char> (in)), std::istreambuf_iterator<char>()));
}

constexpr auto sample_rate = 44100;
//------------------------------------------------------------------------------
class test_playhead : public juce::AudioPlayHead {
public:
  void reset()
  {
    pos_info.bpm                = 120;
    pos_info.timeSigNumerator   = 4;
    pos_info.timeSigDenominator = 4;
    /** The current play position, in samples from the start of the timeline.
     */
    pos_info.timeInSamples = 0;
    /** The current play position, in seconds from the start of the timeline.
     */
    pos_info.timeInSeconds  = 0;
    pos_info.editOriginTime = 0;
    /** The current play position, in units of quarter-notes. */
    pos_info.ppqPosition               = 0;
    pos_info.ppqPositionOfLastBarStart = 0;
    FrameRateType frameRate;
    pos_info.isPlaying    = true;
    pos_info.isRecording  = false;
    pos_info.ppqLoopStart = 0;
    pos_info.ppqLoopEnd   = 0;
    pos_info.isLooping    = false;
  }

  void add_samples (uint n)
  {
    pos_info.timeInSamples += n;
    pos_info.timeInSeconds = ((double) pos_info.timeInSamples) / sample_rate;
    pos_info.ppqPosition   = (pos_info.timeInSeconds / 60.) * pos_info.bpm * 4.;
  }

  juce::AudioPlayHead::CurrentPositionInfo pos_info {};
  bool getCurrentPosition (CurrentPositionInfo& result) override
  {
    result = pos_info;
    return true;
  };
};
//------------------------------------------------------------------------------
class smoke_test : public ::testing::Test {
public:
  smoke_test() {}

  void SetUp()
  {
    float io_cfg = ((float) 1) / ((float) 0xffff);

    gui_emu.emplace();
    plugin.reset (createPluginFilter());
    plugin->enableAllBuses();
    plugin->setPlayHead (&playhead);
    playhead.reset();
    audio.setSize (16, 128);
    plugin->prepareToPlay (sample_rate, 1024);

    for (juce::AudioProcessorParameter* param : plugin->getParameters()) {
      juce::String name = param->getName (64);
      // dbg ("parameter: %s\n", name.toRawUTF8());
      if (name == "fx_type_01") {
        fx_type_param = param;
      }
      else if (name == "in_selection_01") {
        ins_param = param;
        param->setValueNotifyingHost (io_cfg);
        dbg (
          "Ins set to: %s\n",
          param->getText (param->getValue(), 0).toRawUTF8());
      }
      else if (name == "out_selection_01") {
        outs_param = param;
        param->setValueNotifyingHost (io_cfg);
        dbg (
          "Outs set to: %s\n",
          param->getText (param->getValue(), 0).toRawUTF8());
      }
    }
    ASSERT_TRUE (fx_type_param != nullptr);
    ASSERT_TRUE (ins_param != nullptr);
    ASSERT_TRUE (outs_param != nullptr);
  }

  void TearDown()
  {
    plugin->releaseResources();
    plugin.reset();
    gui_emu.reset();
  }

  void fill_audio_buffers_with_noise (float max_level)
  {
    auto l = audio.getWritePointer (0);
    auto r = audio.getWritePointer (1);
    for (uint i = 0; i < audio.getNumSamples(); ++i, ++l, ++r) {
      *l = noise()[0] * max_level;
      *r = *l;
    }
  }

  std::optional<juce::ScopedJuceInitialiser_GUI> gui_emu;
  std::unique_ptr<juce::AudioProcessor>          plugin;
  juce::AudioBuffer<float>                       audio;
  juce::AudioProcessorParameter*                 fx_type_param = nullptr;
  juce::AudioProcessorParameter*                 ins_param     = nullptr;
  juce::AudioProcessorParameter*                 outs_param    = nullptr;
  juce::MidiBuffer                               midi;
  test_playhead                                  playhead;
  white_noise_generator                          noise;
};
//------------------------------------------------------------------------------
#if 0
//------------------------------------------------------------------------------
TEST_F (smoke_test, test_some_preset)
{
  // TODO: probably convert this in a CLI utility...
  std::string patch = get_file_contents ("SOMEPRESET.vstpreset");

  plugin->setStateInformation ((void const*) patch.c_str(), patch.size());

  playhead.reset();

  auto cycles = (sample_rate * 500) / audio.getNumSamples();

  for (uint j = 0; j < cycles; ++j) {
    fill_audio_buffers_with_noise (0.1f);
    plugin->processBlock (audio, midi);

    if ((j & 7) == 0) {
      playhead.add_samples (audio.getNumSamples() / 2);
    }
    else {
      playhead.add_samples (audio.getNumSamples());
    }

    auto l = audio.getWritePointer (0);
    for (uint i = 0; i < audio.getNumSamples(); ++i, ++l) {
      ASSERT_LT (*l, 10.); // spikes, feedback loops, etc.
    }
  }
}
#endif
//------------------------------------------------------------------------------
#define ENGINE_PROFILING 0
#if ENGINE_PROFILING
//------------------------------------------------------------------------------
TEST_F (smoke_test, with_no_fx_on_ch1)
{
  playhead.reset();
  fx_type_param->setValueNotifyingHost (0);
  dbg (
    "Type set to: %s\n",
    fx_type_param->getText (fx_type_param->getValue(), 0).toRawUTF8());

  for (uint i = 0; i < 1000000; ++i) {
    plugin->processBlock (audio, midi);
  }
}
//------------------------------------------------------------------------------
#else
TEST_F (smoke_test, run_all_fx_on_ch1)
{
  // this is mostly to be compiled with a memory sanitizer to be able to
  // verify
  float interval = 1. / fx_type_param->getNumSteps();
  for (float v = 0.; v <= 1.; v += interval) {
    playhead.reset();
    fx_type_param->setValueNotifyingHost (v);
    dbg (
      "Type set to: %s\n",
      fx_type_param->getText (fx_type_param->getValue(), 0).toRawUTF8());
    plugin->processBlock (audio, midi);
  }
}
//------------------------------------------------------------------------------
TEST_F (smoke_test, signal_not_blowing_up)
{
#if 1
  // This is actually testing the whole reserved range.
  // float interval = 1. / 50;
  float interval = 1. / fx_type_param->getNumSteps();
#else
  // For deeply testing some FX
  float interval = 1.;
#endif

  for (float v = 0.; v <= 1.; v += interval) {
    playhead.reset();
    fx_type_param->setValueNotifyingHost (v);
    dbg (
      "Type set to: %s\n",
      fx_type_param->getText (fx_type_param->getValue(), 0).toRawUTF8());

#if 1
    auto cycles = (sample_rate * 5) / audio.getNumSamples();
#else
    // For deeply testing some FX
    auto cycles = (sample_rate * 500) / audio.getNumSamples();
#endif
    for (uint j = 0; j < cycles; ++j) {
      fill_audio_buffers_with_noise (0.1f);
      plugin->processBlock (audio, midi);
      playhead.add_samples (audio.getNumSamples());

      auto l = audio.getWritePointer (0);
      for (uint i = 0; i < audio.getNumSamples(); ++i, ++l) {
        ASSERT_LT (*l, 10.); // spikes, feedback loops, etc.
      }
    }
  }
}
#endif
//------------------------------------------------------------------------------
} // namespace artv
