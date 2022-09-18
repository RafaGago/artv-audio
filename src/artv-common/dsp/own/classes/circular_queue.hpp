#pragma once

#include <cassert>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
//------------------------------------------------------------------------------
// add operators, iterators and stuff if/when needed
template <class T>
class static_pow2_circular_queue {
public:
  //----------------------------------------------------------------------------
  void reset (xspan<T> mem)
  {
    assert (is_pow2 (mem.size()));
    xspan_memset (mem, 0);
    _mem  = mem.data();
    _mask = mem.size() - 1;
    _head = _tail = 0;
  }
  //----------------------------------------------------------------------------
  T pop()
  {
    assert (size());
    auto ret = _mem[_tail & _mask];
    ++_tail;
    return ret;
  }
  //----------------------------------------------------------------------------
  void push (T const& v)
  {
    // the interface could return a bool or throw...
    assert (size() < capacity());
    _mem[_head & _mask] = v;
    ++_head;
  }
  //----------------------------------------------------------------------------
  void push (xspan<const T> vs)
  {
    assert ((vs.size() + size()) <= capacity());
    for (auto v : vs) {
      push (v);
    }
  }
  //----------------------------------------------------------------------------
  uint size() const { return _head - _tail; }
  //----------------------------------------------------------------------------
  uint capacity() const { return _mask + 1; }
  //----------------------------------------------------------------------------
private:
  T*   _mem;
  uint _mask;
  uint _head;
  uint _tail;
};
//------------------------------------------------------------------------------
} // namespace artv
