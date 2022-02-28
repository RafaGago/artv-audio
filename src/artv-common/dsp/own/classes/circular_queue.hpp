#pragma

#include <cassert>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {
//------------------------------------------------------------------------------
// add operators, iterators and stuff if/when needed
template <class T>
class static_pow2_circular_queue {
public:
  //----------------------------------------------------------------------------
  void reset (crange<T> mem)
  {
    assert (is_pow2 (mem.size()));
    crange_memset (mem, 0);
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
    assert (size() <= _mask);
    _mem[_head & _mask] = v;
    ++_head;
  }
  //----------------------------------------------------------------------------
  uint size() const { return _head - _tail; }
  //----------------------------------------------------------------------------
private:
  T*   _mem;
  uint _mask;
  uint _head;
  uint _tail;
};
//------------------------------------------------------------------------------
} // namespace artv
