#pragma once

#include <math.h>
#include <stdint.h>

#include "artv-common/misc/short_ints.hpp"

namespace artv {

struct plugin_play_state {
  double bpm;
  u64    sample_position;
  double quarter_note_position;
  double time_sec_position; // seconds
  bool   is_playing;
  bool   is_valid;
};

class plugin_context {
public:
  virtual ~plugin_context() {};

  virtual uint              get_sample_rate() const               = 0;
  virtual uint              get_max_block_samples() const         = 0;
  virtual plugin_play_state get_play_state() const                = 0;
  virtual void              set_delay_compensation (uint samples) = 0;

  float get_samples_per_beat() const
  {
    return get_play_state().bpm * (1. / 60.) * ((float) get_sample_rate());
  }
};

} // namespace artv
