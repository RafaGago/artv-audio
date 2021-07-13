#pragma once

#include <algorithm>
#include <cstring>
#include <type_traits>
#include <vector>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
// -----------------------------------------------------------------------------
// This class is a helper for doing delay compensation on non-owned buffers,
// buffers that for some reason can't have its size changed, so the only option
// is to save reminders and do memory move operations on them instead of having
// extra space.
// -----------------------------------------------------------------------------
template <class T, uint N_channels = 2, class Index = u32>
class delay_compensation_buffers {
public:
  // ---------------------------------------------------------------------------
  using value_type                 = T;
  static constexpr uint n_channels = N_channels;
  using index_type                 = Index;
  // ---------------------------------------------------------------------------
  template <class U>
  void reset (crange<const U> delays)
  {
    static_assert (std::is_integral_v<U> && std::is_unsigned_v<U>, "");
    // Accepting sample losses when changing delay compensation
    _loc.resize (delays.size() * n_channels);

    uint offset  = 0;
    uint biggest = 0;
    for (uint i = 0; i < delays.size(); ++i) {
      auto delay = delays[i];
      biggest    = std::max<uint> (biggest, delay);
      for (uint j = 0; j < n_channels; ++j) {
        _loc[(i * n_channels) + j].size   = delay;
        _loc[(i * n_channels) + j].offset = offset;
        offset += delay;
      }
    }
    _mem.clear();
    _mem.resize (offset + biggest);
    _tmp_buff = make_crange (&_mem[offset], biggest);
  }
  // ---------------------------------------------------------------------------
  uint size() const noexcept { return _loc.size(); }
  // ---------------------------------------------------------------------------
  void compensate (uint buffer_idx, crange<T> src)
  {
    assert (buffer_idx < size());
    assert (src.size());

    crange<T> dly = get_buffer (buffer_idx);

    if (likely (src.size() >= dly.size())) {
      uint bytes = dly.size() * sizeof (T);
      // save the tail of src samples that are to be delayed on the tmp buffer
      memcpy (&_tmp_buff[0], &src[src.size() - dly.size()], bytes);
      // shift backards the src samples
      memmove (&src[dly.size()], &src[0], (src.size() * sizeof (T)) - bytes);
      // place the delayed samples from the previous round on top
      memcpy (&src[0], &dly[0], bytes);
      // save the src samples that are to be delayed for the next round
      memcpy (&dly[0], &_tmp_buff[0], bytes);
    }
    else {
      uint bytes = src.size() * sizeof (T);
      // save the head of the delayed samples
      memcpy (&_tmp_buff[0], &dly[0], bytes);
      // shift forwards the delayed samples
      memmove (&dly[0], &dly[src.size()], (dly.size() * sizeof (T)) - bytes);
      // move the new samples to the tail of the delayed ones
      memcpy (&dly[dly.size() - src.size()], &src[0], bytes);
      // copy the saved head of the delayed samples to the current sample set
      memcpy (&src[0], &_tmp_buff[0], bytes);
    }
  }
  // ---------------------------------------------------------------------------
private:
  // ---------------------------------------------------------------------------
  crange<T> get_buffer (uint idx)
  {
    return {&_mem[_loc[idx].offset], _loc[idx].size};
  }
  // ---------------------------------------------------------------------------
  struct buffer_location {
    index_type size   = 0;
    index_type offset = 0;
  };
  std::vector<buffer_location> _loc;
  std::vector<T>               _mem;
  crange<T>                    _tmp_buff;
};
// -----------------------------------------------------------------------------
// This class is a helper for doing delay compensation an owned buffers. The
// owned buffer is ironically not owned by the class itself, the class just
// manages them.
//
// The usage is basically calling:
// -get_write_buffer, process and modify there
// -get_read_buffer: to pass the delay compensated buffer somewhere.
// -iteration_end: To prepare the buffer head for the next call
//
// The "delay" is a fixed paramter, but not owned by the class.
// -----------------------------------------------------------------------------
template <class T>
class owned_delay_compensation_buffer {
public:
  //----------------------------------------------------------------------------
  void reset (crange<T> buffer, uint max_delay)
  {
    _buff      = buffer;
    _max_delay = max_delay;
    memset (_buff.data(), 0, sizeof _buff[0] * _buff.size());
  }
  //----------------------------------------------------------------------------
  void iteration_end (uint delay, uint total_bytes_written)
  {
    assert (total_bytes_written <= _buff.size());
    memcpy (
      &_buff[_max_delay - delay],
      &_buff[_max_delay - delay + total_bytes_written],
      delay * sizeof _buff[0]);
  }
  //----------------------------------------------------------------------------
  crange<T> get_write_buffer()
  {
    return {&_buff[_max_delay], _buff.size() - _max_delay};
  }
  //----------------------------------------------------------------------------
  crange<T> get_read_buffer (uint delay)
  {
    return {&_buff[_max_delay - delay], _buff.size() - _max_delay};
  }
  //----------------------------------------------------------------------------
private:
  crange<T> _buff;
  uint      _max_delay;
  // ---------------------------------------------------------------------------
};

} // namespace artv
