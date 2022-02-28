#pragma once

#include <array>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

#include "artv-common/dsp/own/classes/mufft.hpp"
#include "artv-common/dsp/own/classes/pffft.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
//------------------------------------------------------------------------------
namespace detail {
//------------------------------------------------------------------------------
template <class T>
using fft_impl = std::conditional_t<
  std::is_same_v<T, float>,
  mufft::fft<float>,
  pffft::fft<double>>;
//------------------------------------------------------------------------------
template <class FFT>
static constexpr bool fft_requires_buffer = !FFT::io_can_alias
  || FFT::reorder_requires_buffer || FFT::transform_requires_buffer;
//------------------------------------------------------------------------------
template <class T, uint Align, bool has_buffer, bool buffer_is_external>
class fft_members {
public:
  void resize (uint v) { _size = v; }
  uint buffer_size() const { return _size; }

private:
  uint _size;
};
//------------------------------------------------------------------------------
template <class T, uint Align>
class fft_members<T, Align, true, false> {
public:
  void      resize (uint v) { _work.resize (v); }
  uint      buffer_size() const { return _work.size(); }
  crange<T> get_buffer() { return _work; }

private:
  std::vector<T, overaligned_allocator<T, Align>> _work;
};
//------------------------------------------------------------------------------
template <class T, uint Align>
class fft_members<T, Align, true, true> {
public:
  void      resize (uint v) {}
  uint      buffer_size() const { return _work.size(); }
  crange<T> get_buffer() { return _work; }
  void      set_buffer (crange<T> v, uint iosize)
  {
    assert (v.size() >= iosize);
    _work = make_crange (v.data(), iosize);
  }

private:
  crange<T> _work;
};
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
template <class T, bool Use_external_buffer = false>
class fft {
private:
  static_assert (std::is_floating_point_v<T>);
  using impl = detail::fft_impl<T>;
  //----------------------------------------------------------------------------
public:
  using value_type                   = T;
  static constexpr uint io_alignment = impl::io_alignment;
  using allocator = overaligned_allocator<value_type, io_alignment>;
  //----------------------------------------------------------------------------
  bool reset (
    uint               blocksize,
    bool               complex = true,
    crange<value_type> extbuff = {})
  {
    uint iosize = blocksize * (complex ? 2 : 1);
    if (!_impl.reset (blocksize, complex)) {
      return false;
    }
    _work.resize (iosize);
    if constexpr (Use_external_buffer) {
      _work.set_buffer (extbuff, iosize);
    }
    _complex = complex;
    return true;
  }
  //----------------------------------------------------------------------------
  void forward (crange<value_type> out, const crange<value_type> in)
  {
    assert (out.size() >= _work.buffer_size());
    assert (in.size() >= _work.buffer_size());
    if constexpr (impl::transform_requires_buffer) {
      _impl.forward (out, in, _work.get_buffer());
    }
    else {
      _impl.forward (out, in);
    }
  }
  //----------------------------------------------------------------------------
  void forward (crange<value_type> io)
  {
    if constexpr (impl::io_can_alias) {
      forward (io, io);
    }
    else {
      crange_copy (_work.get_buffer(), io);
      forward (io, _work.get_buffer());
    }
  }
  //----------------------------------------------------------------------------
  void forward_ordered (crange<value_type> out, const crange<value_type> in)
  {
    assert (out.size() >= _work.buffer_size());
    assert (in.size() >= _work.buffer_size());
    if constexpr (impl::transform_requires_buffer) {
      _impl.forward_ordered (out, in, _work.get_buffer());
    }
    else {
      _impl.forward_ordered (out, in);
    }
  }
  //----------------------------------------------------------------------------
  void forward_ordered (crange<value_type> io)
  {
    if constexpr (impl::io_can_alias) {
      forward_ordered (io, io);
    }
    else {
      crange_copy (_work.get_buffer(), io);
      forward_ordered (io, _work.get_buffer());
    }
  }
  //----------------------------------------------------------------------------
  void reorder_after_forward (crange<value_type> out)
  {
    if constexpr (impl::reorder_requires_buffer) {
      _impl.reorder_after_forward (out, _work.get_buffer());
    }
    else {
      _impl.reorder_after_forward (out);
    }
  }
  //----------------------------------------------------------------------------
  void backward (crange<value_type> out, const crange<value_type> in)
  {
    assert (out.size() >= _work.buffer_size());
    assert (in.size() >= _work.buffer_size());
    if constexpr (impl::transform_requires_buffer) {
      _impl.backward (out, in, _work.get_buffer());
    }
    else {
      _impl.backward (out, in);
    }
  }
  //----------------------------------------------------------------------------
  void backward (crange<value_type> io)
  {
    if constexpr (impl::io_can_alias) {
      backward (io, io);
    }
    else {
      crange_copy (_work.get_buffer(), io);
      backward (io, _work.get_buffer());
    }
  }
  //----------------------------------------------------------------------------
  void backward_ordered (crange<value_type> out, const crange<value_type> in)
  {
    assert (out.size() >= _work.buffer_size());
    assert (in.size() >= _work.buffer_size());
    if constexpr (impl::transform_requires_buffer) {
      _impl.backward_ordered (out, in, _work.get_buffer());
    }
    else {
      _impl.backward_ordered (out, in);
    }
  }
  //----------------------------------------------------------------------------
  void backward_ordered (crange<value_type> io)
  {
    if constexpr (impl::io_can_alias) {
      backward_ordered (io, io);
    }
    else {
      crange_copy (_work.get_buffer(), io);
      backward_ordered (io, _work.get_buffer());
    }
  }
  //----------------------------------------------------------------------------
  void reorder_before_backward (crange<value_type> out)
  {
    if constexpr (impl::reorder_requires_buffer) {
      _impl.reorder_before_backward (out, _work.get_buffer());
    }
    else {
      _impl.reorder_before_backward (out);
    }
  }
  //----------------------------------------------------------------------------
  void data_rescale (crange<value_type> v)
  {
    _impl.data_rescale (v, size(), is_complex());
  }
  //----------------------------------------------------------------------------
  uint size() const
  {
    return _complex ? _work.buffer_size() / 2 : _work.buffer_size();
  }
  //----------------------------------------------------------------------------
  uint buffer_size() const { return _work.buffer_size(); }
  //----------------------------------------------------------------------------
  bool is_complex() const { return _complex; }
  //----------------------------------------------------------------------------
private:
  static constexpr bool has_work_buffer = detail::fft_requires_buffer<impl>;

  impl _impl;
  detail::
    fft_members<value_type, io_alignment, has_work_buffer, Use_external_buffer>
       _work;
  bool _complex;
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// An ad-hoc class to have multiple FFTs of different sizes initialized, mainly
// to wrap the JSFX fft.
template <class T, bool Complex, size_t... blocksizes>
class initialized_ffts {
public:
  static_assert (std::is_floating_point_v<T>);
  using fft_type                     = fft<T, true>;
  using value_type                   = T;
  static constexpr uint io_alignment = fft_type::io_alignment;
  static constexpr bool is_complex   = Complex;
  //----------------------------------------------------------------------------
  initialized_ffts()
  {
    _workbuff.resize (biggest_blocksize() * (is_complex ? 2 : 1));
    mp_foreach_idx (blocksize_list {}, [=] (auto index, auto bsz) {
      _ffts[index.value].reset (bsz.value, is_complex, _workbuff);
    });
  }
  //----------------------------------------------------------------------------
  ~initialized_ffts()                   = default;
  initialized_ffts (initialized_ffts&&) = default;
  initialized_ffts& operator= (initialized_ffts&&) = default;
  initialized_ffts (initialized_ffts const&)       = delete;
  initialized_ffts& operator= (initialized_ffts const&) = delete;
  //----------------------------------------------------------------------------
  fft_type* get_fft (uint blocksize)
  {
    auto idx = get_fft_idx (blocksize);
    return idx >= 0 ? &_ffts[idx] : nullptr;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static int get_fft_idx (uint blocksize)
  {
    int ret = -1;
    mp_foreach_idx (blocksize_list {}, [=, &ret] (auto index, auto bsz) {
      if (bsz.value == blocksize) {
        ret = index.value;
      }
    });
    return ret;
  }
  //----------------------------------------------------------------------------
  static uint biggest_blocksize()
  {
    uint maxbsz = 0;
    mp11::mp_for_each<blocksize_list> ([&] (auto bsz) {
      maxbsz = (bsz.value > maxbsz) ? bsz.value : maxbsz;
    });
    return maxbsz;
  }
  //----------------------------------------------------------------------------
  using blocksize_list = mp11::mp_list_c<int, blocksizes...>;

  std::array<fft_type, sizeof...(blocksizes)>           _ffts;
  std::vector<value_type, typename fft_type::allocator> _workbuff;
};
//------------------------------------------------------------------------------

} // namespace artv
