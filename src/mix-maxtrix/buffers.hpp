#pragma once

#include <algorithm>
#include <limits>
#include <type_traits>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
// -----------------------------------------------------------------------------
// This class gets really confusing, partially because a juce::AudioBuffer is
// actually a collection of buffers. E.g. on an 8 in plugin it has 16 buffers.
template <class T>
struct buffers {
public:
  using value_type = T;
  static_assert (std::is_floating_point<T>::value, "");
  //----------------------------------------------------------------------------
  buffers (size_t channel_count) : _internal {64, (int) channel_count}
  {
    _channel_count = channel_count;
    _buffers[1]    = &_internal;
  }
  //----------------------------------------------------------------------------
  void on_block (juce::AudioBuffer<T>& block)
  {
    assert (block.getNumChannels() == _channel_count);
    _internal.setSize (
      _channel_count, block.getNumSamples(), false, false, true);
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
  void swap (int bus1, int bus2)
  {
    using simd                  = juce::dsp::SIMDRegister<value_type>;
    constexpr size_t simd_elems = simd::SIMDNumElements;

    assert (bus1 != 0);
    assert (bus2 != 0);

    if (unlikely (bus1 == bus2)) {
      return;
    }

    auto simd_swap = [] (std::array<value_type*, 2> v, uint blocks) {
      value_type* end = v[0] + (blocks * simd_elems);
      while (v[0] < end) {
        auto reg1 = simd::fromRawArray (v[0]);
        auto reg2 = simd::fromRawArray (v[1]);
        reg1.copyToRawArray (v[1]);
        reg2.copyToRawArray (v[0]);
        v[0] += simd_elems;
        v[1] += simd_elems;
      }
    };

    auto remainder_swap = [] (std::array<value_type*, 2> v, uint remainder) {
      value_type* end = v[0] + remainder;
      while (v[0] < end) {
        value_type tmp = *v[0];
        *v[0]          = *v[1];
        *v[1]          = tmp;
        ++v[0];
        ++v[1];
      }
    };

    auto chnls1 = get_write_ptrs<N_bus_chnls> (bus1);
    auto chnls2 = get_write_ptrs<N_bus_chnls> (bus2);

    for (int i = 0; i < N_bus_chnls; ++i) {
      block_divide (
        simd::SIMDRegisterSize,
        make_array (chnls1[i], chnls2[i]),
        sample_count(),
        simd_swap,
        remainder_swap);
    }
  }
  //----------------------------------------------------------------------------
  template <size_t N_bus_chnls = 2, class I> // 2 = stereo
  bool mix (
    int             dst_bus,
    crange<const I> mix_buses,
    bool            foce_sum_with_dst = false)
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
    int             dst_bus,
    crange<const I> mix_buses,
    crange<P*>      summing_procesors,
    bool            foce_sum_with_dst = false)
  {
    static_assert (
      std::is_integral<I>::value && std::is_signed<I>::value
        && std::numeric_limits<I>::max() <= std::numeric_limits<int>::max(),
      "Invalid input array type");
    assert (summing_procesors.size() >= mix_buses.size());

    size_t samples = sample_count();
    uint   mixed   = foce_sum_with_dst ? 1 : 0;
    // check if dst_bus is in the mix sources, it has to go first before writing
    // starts.
    for (uint i = 0; i < mix_buses.size(); ++i) {
      if (dst_bus == mix_buses[i]) {
        auto in     = get_read_ptrs<N_bus_chnls> (dst_bus);
        auto out    = get_write_ptrs<N_bus_chnls> (dst_bus);
        bool do_sum = mixed > 0;
        summing_procesors[i]->process (out, in, samples, do_sum);
        ++mixed;
        break;
      }
    }

    for (uint bus : mix_buses) {
      if (bus != 0 && bus != dst_bus) {
        auto out    = get_write_ptrs<N_bus_chnls> (dst_bus);
        auto in     = get_read_ptrs<N_bus_chnls> (bus);
        bool do_sum = mixed > 0;
        summing_procesors[bus - 1]->process (out, in, samples, do_sum);
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
    using simd                  = juce::dsp::SIMDRegister<value_type>;
    constexpr size_t simd_elems = simd::SIMDNumElements;

    if (unlikely (bus == 0)) {
      return;
    }

    uint first = get_first_audiobuffer_chnl (bus, N_bus_chnls);
    for (int i = 0; i < N_bus_chnls; ++i) {
      if (from[i] == to[i]) {
        _buffers[bus < 0]->applyGain (first + i, 0, sample_count(), from[i]);
        continue;
      }

      auto ptrs = get_write_ptrs<N_bus_chnls> (bus);

      // maybe add interpolation classes?
      auto simd_ramp = [&] (std::array<value_type*, 1> v, uint blocks) {
        value_type* end     = v[0] + (blocks * simd_elems);
        value_type  current = value_type {from[i]};

        value_type increment = value_type {to[i]} - value_type {from[i]};
        increment *= (value_type) simd_elems;
        increment /= (value_type) sample_count();

        while (v[0] < end) {
          auto reg = simd::fromRawArray (v[0]);
          reg *= current;
          reg.copyToRawArray (v[0]);
          v[0] += simd_elems;
          current += increment;
        }
      };

      auto remainder_ramp = [&] (std::array<value_type*, 1> v, uint remainder) {
        value_type* end = v[0] + remainder;
        while (v[0] < end) {
          *v[0] *= to[i];
          ++v[0];
        }
      };

      block_divide (
        simd::SIMDRegisterSize,
        make_array (ptrs[i]),
        sample_count(),
        simd_ramp,
        remainder_ramp);
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
}; // namespace artv

} // namespace artv
