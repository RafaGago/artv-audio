#pragma once

// Copy Boost's
// https://www.boost.org/doc/libs/1_77_0/boost/align/overaligned_allocator.hpp

#include <cstddef>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>

#include "artv-common/misc/bits.hpp"

namespace artv {
//------------------------------------------------------------------------------
template <class T, std::size_t Alignment>
class overaligned_allocator {
public:
  //----------------------------------------------------------------------------
  static constexpr std::size_t alignment = Alignment;
  //----------------------------------------------------------------------------
  using value_type         = T;
  using pointer            = T*;
  using const_pointer      = const T*;
  using void_pointer       = void*;
  using const_void_pointer = const void*;
  using reference          = std::add_lvalue_reference_t<T>;
  using const_reference    = std::add_lvalue_reference_t<const T>;
  using size_type          = std::size_t;
  using difference_type    = std::ptrdiff_t;
  using propagate_on_container_move_assignment = std::true_type;
  using is_always_equal                        = std::true_type;
  //----------------------------------------------------------------------------
  static_assert (alignment >= alignof (T), "Underaligning is not supported");
  static_assert (is_pow2 (alignment), "Aligment is not a power of 2");
  //----------------------------------------------------------------------------
  template <class U>
  struct rebind {
    typedef overaligned_allocator<U, alignment> other;
  };
  //----------------------------------------------------------------------------
  overaligned_allocator() = default;
  //----------------------------------------------------------------------------
  template <class U>
  overaligned_allocator (overaligned_allocator<U, alignment> const&) noexcept
  {}
  //----------------------------------------------------------------------------
  pointer allocate (size_type size, const_void_pointer = 0)
  {
    void* ret;
    if constexpr (alignof (std::max_align_t) >= alignment) {
      ret = operator new (sizeof (T) * size);
    }
    else {
      ret = operator new (sizeof (T) * size, std::align_val_t {alignment});
    }
    return static_cast<T*> (ret);
  }
  //----------------------------------------------------------------------------
  void deallocate (pointer ptr, size_type)
  {
    if constexpr (alignof (std::max_align_t) >= alignment) {
      operator delete (ptr);
    }
    else {
      operator delete (ptr, std::align_val_t {alignment});
    }
  }
  //----------------------------------------------------------------------------
  constexpr size_type max_size() const noexcept
  {
    return std::numeric_limits<size_type>::max() / sizeof (value_type);
  }
  //----------------------------------------------------------------------------
  template <class U, class... Args>
  void construct (U* ptr, Args&&... args)
  {
    ::new ((void*) ptr) U (std::forward<Args> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class U>
  void destroy (U* ptr)
  {
    (void) ptr;
    ptr->~U();
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <class T, class U, std::size_t Alignment>
inline bool operator== (
  const overaligned_allocator<T, Alignment>&,
  const overaligned_allocator<U, Alignment>&) noexcept
{
  return true;
}
//------------------------------------------------------------------------------
template <class T, class U, std::size_t Alignment>
inline bool operator!= (
  const overaligned_allocator<T, Alignment>&,
  const overaligned_allocator<U, Alignment>&) noexcept
{
  return false;
}

//------------------------------------------------------------------------------
} // namespace artv
