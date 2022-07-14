#pragma once

#include <cassert>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"

namespace artv {

//------------------------------------------------------------------------------
class juce_plugin_context : public plugin_context {
public:
  ~juce_plugin_context() {};
  //------------------------------------------------------------------------------
  void reset (juce::AudioProcessor& proc, uint samplerate, uint max_block_size)
  {
    _samplerate     = samplerate;
    _max_block_size = max_block_size;
    _processor      = &proc;
  }
  //----------------------------------------------------------------------------
  uint get_sample_rate() const override { return _samplerate; }
  //----------------------------------------------------------------------------
  uint get_max_block_samples() const override { return _max_block_size; }
  //----------------------------------------------------------------------------
  virtual plugin_play_state get_play_state() const override
  {
    auto pos = _processor->getPlayHead()->getPosition();
    if (pos) {
      auto bpm    = pos->getBpm();
      auto t_spls = pos->getTimeInSamples();
      auto ppq    = pos->getPpqPosition();
      auto sec    = pos->getTimeInSeconds();
      return {
        .bpm                   = bpm ? *bpm : 0.,
        .sample_position       = t_spls ? (u64) *t_spls : 0,
        .quarter_note_position = ppq ? *ppq : 0,
        .time_sec_position     = sec ? *sec : 0,
        .is_playing            = pos->getIsPlaying(),
        .is_valid              = bpm && t_spls && ppq && sec,
      };
    };
    return {};
  }
  //----------------------------------------------------------------------------
  void set_delay_compensation (uint samples) override
  {
    _processor->setLatencySamples (samples);
  }
  //------------------------------------------------------------------------------
private:
  juce::AudioProcessor* _processor      = nullptr;
  uint                  _samplerate     = 0;
  uint                  _max_block_size = 0;
};
//------------------------------------------------------------------------------
} // namespace artv
