// own old invention of (a dynamic extent) std::span, (previously named xspan)
// prepared from on C++20 probably inherit an std::span<T, std::dynamic_extent>
#pragma once

#include <array>
#include <vector>

#include "artv-common/misc/short_ints.hpp"

namespace artv {
//------------------------------------------------------------------------------
// TODO: this won't match a std::span, it is just the old own version I was
// using.
template <class T>
class xspan {
private:
  //----------------------------------------------------------------------------
  template <class U>
  static constexpr bool same_or_non_const_to_const
    = std::is_same_v<T, U> || std::is_same_v<std::remove_const_t<T>, U>;
  //----------------------------------------------------------------------------
  // This one gets "U" just to have a dependant type for SFINAE. It decides
  // based on "T".
  template <class U>
  static constexpr bool      xspan_type_is_const
    = std::is_same_v<U, U>&& std::is_const_v<T>;
  //----------------------------------------------------------------------------
public:
  using value_type      = T;
  using iterator        = value_type*;
  using const_iterator  = value_type const*;
  using pointer         = value_type*;
  using const_pointer   = value_type const*;
  using reference       = value_type&;
  using const_reference = value_type const&;

  constexpr xspan() = default;

  constexpr xspan (value_type* start, size_t size)
  {
    _start = start;
    _size  = start ? size : 0;
  }

  template <class U, std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr xspan (xspan<U> const& other) : xspan {other.data(), other.size()}
  {}

  template <class U, std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr xspan<T>& operator= (xspan<U> const& other)
  {
    _start = other.data();
    _size  = other.size();
    return *this;
  }

  template <
    class U,
    uint N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr xspan (U (&arr)[N]) : xspan {&arr[0], N}
  {}

  template <
    class U,
    uint N,
    std::enable_if_t<xspan_type_is_const<U>>* = nullptr>
  constexpr xspan (U const (&arr)[N]) : xspan {&arr[0], N}
  {}

  template <
    class U,
    uint N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr xspan<T>& operator= (U (&arr)[N])
  {
    _start = &arr[0];
    _size  = N;
    return *this;
  }

  template <
    class U,
    uint N,
    std::enable_if_t<xspan_type_is_const<U>>* = nullptr>
  constexpr xspan<T>& operator= (U const (&arr)[N])
  {
    _start = &arr[0];
    _size  = N;
    return *this;
  }

  template <
    class U,
    size_t N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr xspan (std::array<U, N>& arr) : xspan {arr.data(), arr.size()}
  {}

  template <
    class U,
    size_t N,
    std::enable_if_t<xspan_type_is_const<U>>* = nullptr>
  constexpr xspan (std::array<U, N> const& arr) : xspan {arr.data(), arr.size()}
  {}

  template <class U, size_t N>
  xspan (std::array<U, N>&& arr)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
  }

  template <
    class U,
    size_t N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr xspan<T>& operator= (std::array<U, N>& arr)
  {
    _start = arr.data();
    _size  = N;
    return *this;
  }

  template <
    class U,
    size_t N,
    std::enable_if_t<xspan_type_is_const<U>>* = nullptr>
  constexpr xspan<T>& operator= (std::array<U, N> const& arr)
  {
    _start = arr.data();
    _size  = N;
    return *this;
  }

  template <class U, size_t N>
  xspan<T>& operator= (std::array<U, N>&& arr)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
    return *this;
  }

  // remember that xspan is a non-owning reference that can get
  // easily invalidated...
  template <
    class U,
    class Alloc,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  xspan (std::vector<U, Alloc>& vec) : xspan {vec.data(), vec.size()}
  {}

  template <
    class U,
    class Alloc,
    std::enable_if_t<xspan_type_is_const<U>>* = nullptr>
  xspan (std::vector<U, Alloc> const& vec) : xspan {vec.data(), vec.size()}
  {}

  template <class U, class Alloc>
  xspan (std::vector<U, Alloc>&& vec)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
  }

  template <
    class U,
    class Alloc,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  xspan<T>& operator= (std::vector<U, Alloc>& vec)
  {
    _start = vec.data();
    _size  = vec.size();
    return *this;
  }

  template <
    class U,
    class Alloc,
    std::enable_if_t<xspan_type_is_const<U>>* = nullptr>
  xspan<T>& operator= (std::vector<U, Alloc>& vec)
  {
    _start = vec.data();
    _size  = vec.size();
    return *this;
  }

  template <class U, class Alloc>
  xspan<T>& operator= (std::vector<U, Alloc>&& vec)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
    return *this;
  }

  constexpr value_type& at (size_t idx)
  {
    assert (_start);
    assert (idx < size());
    return _start[idx];
  }

  constexpr value_type const& at (size_t idx) const
  {
    assert (_start);
    assert (idx < size());
    return _start[idx];
  }

  constexpr value_type&       operator[] (size_t idx) { return at (idx); }
  constexpr value_type const& operator[] (size_t idx) const { return at (idx); }

  constexpr explicit operator bool() const { return !empty(); }

  // returns a copy/subrange with "count" elements dropped from the tail.
  constexpr xspan<T> reduced (uint count) const
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._size -= count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }

  // returns a copy/subrange with "count" elements dropped from the head.
  constexpr xspan<T> advanced (uint count) const
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._start += count;
    r._size -= count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  // returns a copy/subrange containing "count" elements starting from the
  // head.
  constexpr xspan<T> get_head (uint count) const
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._size  = count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  // returns a copy/subrange containing "count" elements starting from the tail.
  constexpr xspan<T> get_tail (uint count) const
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._start += r._size - count;
    r._size  = count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  // drops "count" elems from the head and returns the cut subrange
  constexpr xspan<T> cut_head (uint count)
  {
    auto r = get_head (count);
    *this  = advanced (count);
    return r;
  }
  // drops "count" elems from the tail and returns the cut subrange
  constexpr xspan<T> cut_tail (uint count)
  {
    auto r = get_tail (count);
    *this  = reduced (count);
    return r;
  }

  constexpr iterator       begin() { return _start; }
  constexpr const_iterator cbegin() { return _start; }
  constexpr const_iterator begin() const { return _start; }

  constexpr iterator       end() { return _start + _size; }
  constexpr const_iterator cend() { return _start + _size; }
  constexpr const_iterator end() const { return _start + _size; }

  constexpr size_t   size() const { return _size; }
  constexpr size_t   byte_size() const { return _size * sizeof (value_type); }
  constexpr T*       data() { return _start; }
  constexpr T const* data() const { return _start; }

  constexpr bool empty() const { return size() == 0; }

  constexpr T&       first() { return at (0); }
  constexpr T const& first() const { return at (0); }
  constexpr T&       last() { return at (_size - 1); }
  constexpr T const& last() const { return at (_size - 1); }

  constexpr void clear()
  {
    _start = nullptr;
    _size == 0;
  }

  template <class U>
  xspan<U> cast() const
  {
    constexpr auto big   = std::max (sizeof (value_type), sizeof (U));
    constexpr auto small = std::min (sizeof (value_type), sizeof (U));

    static_assert ((big % small) == 0, "sizes are not multiples");

    return {reinterpret_cast<U*> (_start), byte_size() / sizeof (U)};
  }
  // dummy parameter version to avoid calling make_xspan(x).template cast<T>();
  template <class U>
  xspan<U> cast (U) const
  {
    return cast<U>();
  }

  xspan<std::add_const_t<T>> to_const() const { return *this; }

private:
  //----------------------------------------------------------------------------
  T*     _start = nullptr;
  size_t _size  = 0;
};
//------------------------------------------------------------------------------
template <class T>
static constexpr xspan<T> make_xspan (T* mem, size_t count)
{
  return {mem, count};
};

template <class T>
static constexpr xspan<const T> make_xspan (T const* mem, size_t count)
{
  return {mem, count};
};

template <class T>
static constexpr xspan<T> make_xspan (T& mem)
{
  return {&mem, 1};
};

template <class T>
static constexpr xspan<const T> make_xspan (T const& mem)
{
  return {&mem, 1};
};

template <class T>
static constexpr void make_xspan (T&& mem)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T, size_t N>
static constexpr xspan<T> make_xspan (
  T (&arr)[N],
  size_t count      = N,
  size_t offset_idx = 0)
{
  assert (count + offset_idx <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static constexpr xspan<const T> make_xspan (
  T const (&arr)[N],
  size_t count      = N,
  size_t offset_idx = 0)
{
  assert (count + offset_idx <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static constexpr xspan<T> make_xspan (
  std::array<T, N>& arr,
  size_t            count      = N,
  size_t            offset_idx = 0)
{
  assert ((count + offset_idx) <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static constexpr xspan<const T> make_xspan (
  std::array<T, N> const& arr,
  size_t                  count      = N,
  size_t                  offset_idx = 0)
{
  assert ((count + offset_idx) <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static constexpr void make_xspan (
  std::array<T, N>&& arr,
  size_t             count      = N,
  size_t             offset_idx = 0)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T, class Alloc>
static constexpr xspan<T> make_xspan (std::vector<T, Alloc>& vec)
{
  return {vec.data(), vec.size()};
};

template <class T, class Alloc>
static constexpr xspan<const T> make_xspan (std::vector<T, Alloc> const& vec)
{
  return {vec.data(), vec.size()};
};

template <class T, class Alloc>
static constexpr void make_xspan (std::vector<T, Alloc>&& vec)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T, class Alloc>
static constexpr xspan<T> make_xspan (
  std::vector<T, Alloc>& vec,
  size_t                 count,
  size_t                 offset_idx = 0)
{
  assert (count + offset_idx <= vec.size() && "out of bounds");
  return {vec.data() + offset_idx, count};
};

template <class T, class Alloc>
static constexpr xspan<const T> make_xspan (
  std::vector<T, Alloc> const& vec,
  size_t                       count,
  size_t                       offset_idx = 0)
{
  assert (count + offset_idx <= vec.size() && "out of bounds");
  return {vec.data() + offset_idx, count};
};

template <class T, class Alloc>
static constexpr void make_xspan (
  std::vector<T, Alloc>&& vec,
  size_t                  count,
  size_t                  offset_idx = 0)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T>
static constexpr xspan<T> make_xspan (
  xspan<T> range,
  size_t   count,
  size_t   offset_idx = 0)
{
  assert (count + offset_idx <= range.size() && "out of bounds");
  return {&range[offset_idx], count};
};

template <class T>
static constexpr xspan<const T> make_xspan (
  xspan<T> range,
  size_t   count,
  size_t   offset_idx = 0)
{
  assert (count + offset_idx <= range.size() && "out of bounds");
  return {&range[offset_idx], count};
};
//------------------------------------------------------------------------------
template <class T>
static void xspan_memset (xspan<T> range, int value = 0)
{
  memset (range.data(), value, range.size() * sizeof range[0]);
}
//------------------------------------------------------------------------------
template <class T, class U>
static uint xspan_memcpy (xspan<T> dst, xspan<const U> src)
{
  uint size = std::min (dst.size() * sizeof (T), src.size() * sizeof (U));
  memcpy (dst.data(), src.data(), size);
  return size;
}
// just for template deduction to work with const
template <class T, class U>
static uint xspan_memcpy (xspan<T> dst, xspan<U> src)
{
  return xspan_memcpy<T, const U> (dst, src);
}
//------------------------------------------------------------------------------
template <class T>
static uint xspan_copy (xspan<T> dst, const xspan<T> src)
{
  return xspan_memcpy (dst, src) / sizeof (T);
}
//------------------------------------------------------------------------------
} // namespace artv
