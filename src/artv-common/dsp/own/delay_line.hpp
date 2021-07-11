#pragma once

#include <array>
#include <cstring>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// just a single allocation multichannel delay line with random access. Capacity
// is always rounded up to a power of two to avoid divisions when wrapping
// around.
template <class T>
class delay_line {
public:
  //----------------------------------------------------------------------------
  void reset (uint channels, uint max_delay)
  {
    reset_with_extra_tail_samples (channels, max_delay, 0);
  }
  //----------------------------------------------------------------------------
  void get (crange<T> channels_out, uint idx) // idx 0 is previously pushed
  {
    assert (channels_out.size() >= _channels);
    assert (idx < max_delay_samples());

    uint pos = ((_wpos + idx) & _capacity_mask) * _channels;
    memcpy (channels_out.data(), &_samples[pos], sizeof (T) * _channels);
  }
  //----------------------------------------------------------------------------
  T get (uint channel, uint idx)
  {
    assert (channel < _channels);
    assert (idx < max_delay_samples());

    return _samples[(((_wpos + idx) & _capacity_mask) * _channels) + channel];
  }
  //----------------------------------------------------------------------------
  void push (crange<const T> channels_in)
  {
    assert (channels_in.size() >= _channels);
    --_wpos;
    uint pos = (_wpos & _capacity_mask) * _channels;
    memcpy (&_samples[pos], channels_in.data(), sizeof (T) * _channels);
  }
  //----------------------------------------------------------------------------
  uint max_delay_samples() const { return _capacity_mask + 1; }
  //----------------------------------------------------------------------------
  uint n_channels() const { return _channels; }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  T* reset_with_extra_tail_samples (
    uint channels,
    uint max_delay,
    uint tail_samples)
  {
    max_delay = pow2_round_ceil (max_delay);
    _samples.clear();
    _samples.resize ((max_delay + tail_samples) * channels); // 0 fill
    _wpos          = 0;
    _capacity_mask = max_delay - 1;
    _channels      = channels;
    return &_samples[_samples.size() - channels];
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  std::vector<T> _samples;
  uint           _wpos;
  uint           _capacity_mask;
  uint           _channels;
};
//------------------------------------------------------------------------------
// untested, wasn't needed after all...
template <class T>
class allpass_interpolated_delay_line : public delay_line<T> {
public:
  //----------------------------------------------------------------------------
  // shadowing.
  //----------------------------------------------------------------------------
  void reset (uint channels, uint capacity_log2)
  {
    _allpass_states
      = this->reset_with_extra_tail_samples (channels, capacity_log2, 1);
  }
  //----------------------------------------------------------------------------
  T get (uint channel, float idx)
  {
    uint  i_idx = (float) idx;
    float fract = idx - i_idx;
    return get (channel, i_idx, get_fractional_coeff (fract));
  }
  //----------------------------------------------------------------------------
  constexpr double get_fractional_coeff (double fraction) // 0 to 1
  {
    return (1.0 - fraction) / (1.0 + fraction);
  }
  //----------------------------------------------------------------------------
  // just to be able to save divisions in cases the fraction is constant. See
  // "get (uint channel, float idx)"
  T get (uint channel, uint idx, double fractional_coeff)
  {
    assert (channel < this->n_channels());
    assert (idx < (this->max_delay_samples() - 1));

    std::array<T, 2> points = {get (channel, idx), get (channel, idx + 1)};
    // this is built based on 1 push 1 retrieval
    return warped_allpass_interpolate (
      fractional_coeff, _allpass_states[channel], points);
  }
  //----------------------------------------------------------------------------
private:
  // magic from RS-MET's library.
  double warped_allpass_interpolate (
    double           coeff,
    float&           prev,
    std::array<T, 2> points)
  {
    double ret = points[0] + (coeff * points[1]) - (0.999 * coeff * prev);
    prev       = ret;
    return ret;
  }
  //----------------------------------------------------------------------------
  T* _allpass_states;
};

} // namespace artv
