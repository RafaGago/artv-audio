#pragma once

#include <algorithm>
#include <limits>
#include <type_traits>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
// -----------------------------------------------------------------------------
// This class gets really confusing, partially because a juce::AudioBuffer is
// actually a collection of buffers. E.g. on an 8 in plugin it has 16 buffers.
//
// On this class the channel index refers to a null channel. Positive indexes
// get the DAW/JUCE owned buffers. Negative indexes the internally owned ones.
template <class T>
struct buffers {
public:
  using value_type = T;
  static_assert (std::is_floating_point<T>::value, "");
  //----------------------------------------------------------------------------
  // channel_count: the number of mono buses required for storing the DAW passed
  // buffers. An equal number of internal buses will be created.
  //
  // extra_count: An aditional number of internal mono buses to create.
  buffers (size_t channel_count, size_t extra_count = 0)
    : _internal {(int) channel_count + (int) extra_count, 64}
  {
    _channel_count = channel_count;
    _extra_count   = extra_count;
    _buffers[1]    = &_internal;
  }
  //----------------------------------------------------------------------------
  void on_block (juce::AudioBuffer<T>& block)
  {
    assert (block.getNumChannels() == _channel_count);
    _internal.setSize (
      _channel_count + _extra_count, block.getNumSamples(), false, false, true);
    _buffers[0] = &block;
  }
  //----------------------------------------------------------------------------
  size_t sample_count() const noexcept { return _buffers[0]->getNumSamples(); }
  //----------------------------------------------------------------------------
  value_type* get_write_ptr (int chnl)
  {
    return get_buffer (chnl, nullptr, [] (juce::AudioBuffer<T>& b, uint c) {
      return b.getWritePointer (c);
    });
  }
  //----------------------------------------------------------------------------
  value_type const* get_read_ptr (int chnl)
  {
    return get_buffer (chnl, nullptr, [] (juce::AudioBuffer<T>& b, uint c) {
      return b.getReadPointer (c);
    });
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2> // 2 = stereo
  static std::array<int, N_bus_chnls> to_mono_buses (int bus)
  {
    int to_next = bus < 0 ? -1 : 1;
    to_next     = (bus != 0) ? to_next : 0;
    int current = (bus - to_next) * N_bus_chnls;

    std::array<int, N_bus_chnls> ret;
    for (int& v : ret) {
      current += to_next;
      v = current;
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  static int to_mono_bus_first (int bus, int n_bus_channels)
  {
    int to_next = bus < 0 ? -1 : 1;
    to_next     = (bus != 0) ? to_next : 0;
    int current = (bus - to_next) * n_bus_channels;
    return current + to_next;
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2> // 2 = stereo
  std::array<value_type*, N_bus_chnls> get_write_ptrs (int bus)
  {
    std::array<value_type*, N_bus_chnls> r {};
    assert (_buffers[0]->getNumChannels() % N_bus_chnls == 0);

    if (unlikely (bus == 0)) {
      return r;
    }

    uint first = get_first_audiobuffer_chnl (bus, N_bus_chnls);
    for (int i = 0; i < N_bus_chnls; ++i) {
      r[i] = get_audiobuffer (bus).getWritePointer (first + i);
    }
    return r;
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2> // 2 = stereo
  std::array<value_type const*, N_bus_chnls> get_read_ptrs (int bus)
  {
    std::array<value_type const*, N_bus_chnls> r {};
    assert (_buffers[0]->getNumChannels() % N_bus_chnls == 0);

    if (unlikely (bus == 0)) {
      return r;
    }

    uint first = get_first_audiobuffer_chnl (bus, N_bus_chnls);
    for (int i = 0; i < N_bus_chnls; ++i) {
      r[i] = get_audiobuffer (bus).getReadPointer (first + i);
    }
    return r;
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2> // 2 = stereo
  void clear (int bus)
  {
    auto ptrs = get_write_ptrs (bus);
    for (value_type* ptr : ptrs) {
      memset (ptr, 0, sample_count() * sizeof *ptr);
    }
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2> // 2 = stereo
  void swap (int bus1, int bus2)
  {
    auto chnls1 = get_write_ptrs<N_bus_chnls> (bus1);
    auto chnls2 = get_write_ptrs<N_bus_chnls> (bus2);
    auto count  = sample_count();
    for (int i = 0; i < N_bus_chnls; ++i) {
      std::swap_ranges (chnls1[i], chnls1[i] + count, chnls2[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2, class I> // 2 = stereo
  bool mix (
    int            dst_bus,
    xspan<I const> mix_buses,
    bool           foce_sum_with_dst = false)
  {
    static_assert (
      std::is_integral<I>::value && std::is_signed<I>::value
        && std::numeric_limits<I>::max() <= std::numeric_limits<int>::max(),
      "Invalid input array type");

    size_t samples = sample_count();
    uint   mixed   = foce_sum_with_dst ? 1 : 0;

    // check if "dst_bus" is in "mix_buses", so we add the first input we find
    // instead of just memcpy'ng.
    for (uint i = 0; i < mix_buses.size(); ++i) {
      if (dst_bus == mix_buses[i]) {
        ++mixed;
        break;
      }
    }

    for (int bus : mix_buses) {
      if (bus != 0 && bus != dst_bus) {
        auto out = get_write_ptrs<N_bus_chnls> (dst_bus);
        auto in  = get_read_ptrs<N_bus_chnls> (bus);
        if (mixed > 0) {
          for (int c = 0; c < N_bus_chnls; ++c) {
            juce::FloatVectorOperations::add (out[c], in[c], samples);
          }
        }
        else {
          for (int c = 0; c < N_bus_chnls; ++c) {
            memcpy (out[c], in[c], samples * sizeof (T));
          }
        }
        ++mixed;
      }
    }

    if (mixed) {
      return true;
    }
    else {
      auto out = get_write_ptrs<N_bus_chnls> (dst_bus);
      for (int c = 0; c < N_bus_chnls; ++c) {
        memset (out[c], 0, samples * sizeof (T));
      }
      return false;
    }
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2, class I, class P> // 2 = stereo
  bool mix (
    int            dst_bus,
    xspan<I const> mix_buses,
    xspan<P*>      adders,
    bool           dst_has_valid_data          = false,
    bool           adder_idx_is_mix_bus_offset = false)
  {
    static_assert (
      std::is_integral<I>::value && std::is_signed<I>::value
        && std::numeric_limits<I>::max() <= std::numeric_limits<int>::max(),
      "Invalid input array type");
    assert (adders.size() >= mix_buses.size());

    size_t samples = sample_count();
    uint   mixed   = dst_has_valid_data ? 1 : 0;
    // check if dst_bus is in the mix sources, it has to go first before
    for (uint i = 0; i < mix_buses.size(); ++i) {
      if (dst_bus == mix_buses[i]) {
        auto in     = get_read_ptrs<N_bus_chnls> (dst_bus);
        auto out    = get_write_ptrs<N_bus_chnls> (dst_bus);
        bool do_sum = mixed > 0;
        adders[i]->process (out, in, samples, do_sum);
        ++mixed;
        break;
      }
    }
    for (uint i = 0; i < mix_buses.size(); ++i) {
      uint bus = mix_buses[i];
      if (bus != 0 && bus != dst_bus) {
        auto out    = get_write_ptrs<N_bus_chnls> (dst_bus);
        auto in     = get_read_ptrs<N_bus_chnls> (bus);
        bool do_sum = mixed > 0;
        // sends are fixed-position, can even have negative indexes. Inputs are
        // always positive.
        uint idx = adder_idx_is_mix_bus_offset ? i : bus - 1;
        assert (bus > 0 || adder_idx_is_mix_bus_offset);
        adders[idx]->process (out, in, samples, do_sum);
        ++mixed;
      }
    }

    if (mixed) {
      return true;
    }
    else {
      auto out = get_write_ptrs<N_bus_chnls> (dst_bus);
      for (int c = 0; c < N_bus_chnls; ++c) {
        memset (out[c], 0, samples * sizeof (T));
      }
      return false;
    }
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2> // 2 = stereo
  void apply_gains (int bus, std::array<float, N_bus_chnls> gains)
  {
    if (unlikely (bus == 0)) {
      return;
    }

    uint first = get_first_audiobuffer_chnl (bus, N_bus_chnls);
    for (int i = 0; i < N_bus_chnls; ++i) {
      get_audiobuffer (bus).applyGain (first + i, 0, sample_count(), gains[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2> // 2 = stereo
  void apply_gains_ramp (
    int                            bus,
    std::array<float, N_bus_chnls> from,
    std::array<float, N_bus_chnls> to)
  {
    if (unlikely (bus == 0)) {
      return;
    }
    uint first = get_first_audiobuffer_chnl (bus, N_bus_chnls);
    for (int i = 0; i < N_bus_chnls; ++i) {
      if (from[i] != to[i]) {
        _buffers[bus < 0]->applyGainRamp (
          first + i, 0, sample_count(), from[i], to[i]);
      }
      else {
        _buffers[bus < 0]->applyGain (first + i, 0, sample_count(), from[i]);
      }
    }
  }
  //----------------------------------------------------------------------------
  juce::AudioBuffer<T>& get_audiobuffer (int chnl)
  {
    assert (_buffers[0]);
    return *_buffers[chnl < 0];
  }
  //----------------------------------------------------------------------------
  uint get_first_audiobuffer_chnl (int bus, size_t n_bus_chnls)
  {
    bus = bus > 0 ? bus : -bus;
    return (bus - 1) * n_bus_chnls;
  }

private:
  //----------------------------------------------------------------------------
  template <class F, class U>
  auto get_buffer (int chnl, U onfail, F func)
  {
    using ret_type = decltype (func (*_buffers[0], 0));
    if (unlikely (chnl == 0)) {
      return ret_type {onfail};
    }
    return func (get_audiobuffer (chnl), get_first_audiobuffer_chnl (chnl, 1));
  }
  //----------------------------------------------------------------------------
  template <class F, class U>
  void get_buffer (int chnl, F func)
  {
    if (unlikely (chnl == 0)) {
      return;
    }
    func (get_audiobuffer (chnl), get_first_audiobuffer_chnl (chnl, 1));
  }
  //----------------------------------------------------------------------------
  // _buffers[1] always points to "_internal". Reminder: a "juce::AudioBuffer"
  // should actually be called "juce::AudioBuffers", as it contains as much
  // buffers as channels the plugin has.
  std::array<juce::AudioBuffer<T>*, 2> _buffers;
  juce::AudioBuffer<T>                 _internal;
  uint                                 _channel_count = 0;
  uint                                 _extra_count   = 0;
};

} // namespace artv
