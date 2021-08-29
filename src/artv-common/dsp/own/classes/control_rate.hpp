#pragma once

#include <cassert>
#include <cstring>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {
//------------------------------------------------------------------------------
class control_rate {
public:
  //----------------------------------------------------------------------------
  constexpr void reset (uint period)
  {
    _pos    = -1;
    _next   = 0;
    _period = period;
  }
  //----------------------------------------------------------------------------
  // returns at which sample of the block the control rate had to act. Negative
  // otherwise.
  constexpr int tick (uint n_samples)
  {
    // this algorithm is based in that the mask is always slower or equal than
    // the block size.
    assert (n_samples <= _period);

    uint pos     = _pos + n_samples;
    _pos         = pos;
    int  expired = (int) (pos - _next);
    bool is_ctrl = expired >= 0;
    _next += is_ctrl ? _period : 0;
    uint offset = n_samples - expired - 1; // always positive
    return is_ctrl ? offset : -1;
  }
  //----------------------------------------------------------------------------
  constexpr uint get_period() const { return _period; }
  //----------------------------------------------------------------------------
private:
  uint _pos;
  uint _next;
  uint _period;
};
//------------------------------------------------------------------------------
} // namespace artv
