#pragma once

#include <cassert>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/plugin_context.hpp"

namespace artv {

class juce_plugin_context : public plugin_context {
public:
  ~juce_plugin_context() {};

  void reset (juce::AudioProcessor& proc, uint samplerate, uint max_block_size)
  {
    _samplerate     = samplerate;
    _max_block_size = max_block_size;
    _processor      = &proc;
  }

  uint get_sample_rate() const override { return _samplerate; }

  uint get_max_block_samples() const override { return _max_block_size; }

  virtual plugin_play_state get_play_state() const override
  {
    juce::AudioPlayHead::CurrentPositionInfo inf;
    if (_processor->getPlayHead()->getCurrentPosition (inf)) {
      return {
        .bpm                   = inf.bpm,
        .sample_position       = (u64) inf.timeInSamples,
        .quarter_note_position = inf.ppqPosition,
        .time_sec_position     = inf.timeInSeconds,
        .is_playing            = inf.isPlaying,
        .is_valid              = true,
      };
    };
    return {};
  }

  void set_delay_compensation (uint samples) override
  {
    _processor->setLatencySamples (samples);
  }

private:
  juce::AudioProcessor* _processor      = nullptr;
  uint                  _samplerate     = 0;
  uint                  _max_block_size = 0;
};

} // namespace artv
