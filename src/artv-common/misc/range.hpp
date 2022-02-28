#pragma once

#include <array>
#include <vector>

#include "artv-common/misc/short_ints.hpp"

namespace artv {
//------------------------------------------------------------------------------
template <class T>
class contiguous_range // std might have something like this. Check.
{
private:
  //----------------------------------------------------------------------------
  template <class U>
  static constexpr bool same_or_non_const_to_const
    = std::is_same_v<T, U> || std::is_same_v<std::remove_const_t<T>, U>;
  //----------------------------------------------------------------------------
  // This one gets "U" just to have a dependant type for SFINAE. It decides
  // based on "T".
  template <class U>
  static constexpr bool      crange_type_is_const
    = std::is_same_v<U, U>&& std::is_const_v<T>;
  //----------------------------------------------------------------------------
public:
  using value_type     = T;
  using iterator       = value_type*;
  using const_iterator = value_type const*;

  constexpr contiguous_range() = default;

  constexpr contiguous_range (value_type* start, size_t size)
  {
    _start = start;
    _size  = start ? size : 0;
  }

  template <class U, std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr contiguous_range (contiguous_range<U> const& other)
    : contiguous_range {other.data(), other.size()}
  {}

  template <class U, std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr contiguous_range<T>& operator= (contiguous_range<U> const& other)
  {
    _start = other.data();
    _size  = other.size();
    return *this;
  }

  template <
    class U,
    uint N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr contiguous_range (U (&arr)[N]) : contiguous_range {&arr[0], N}
  {}

  template <
    class U,
    uint N,
    std::enable_if_t<crange_type_is_const<U>>* = nullptr>
  constexpr contiguous_range (U const (&arr)[N]) : contiguous_range {&arr[0], N}
  {}

  template <
    class U,
    uint N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr contiguous_range<T>& operator= (U (&arr)[N])
  {
    _start = &arr[0];
    _size  = N;
    return *this;
  }

  template <
    class U,
    uint N,
    std::enable_if_t<crange_type_is_const<U>>* = nullptr>
  constexpr contiguous_range<T>& operator= (U const (&arr)[N])
  {
    _start = &arr[0];
    _size  = N;
    return *this;
  }

  template <
    class U,
    size_t N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr contiguous_range (std::array<U, N>& arr)
    : contiguous_range {arr.data(), arr.size()}
  {}

  template <
    class U,
    size_t N,
    std::enable_if_t<crange_type_is_const<U>>* = nullptr>
  constexpr contiguous_range (std::array<U, N> const& arr)
    : contiguous_range {arr.data(), arr.size()}
  {}

  template <class U, size_t N>
  contiguous_range (std::array<U, N>&& arr)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
  }

  template <
    class U,
    size_t N,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  constexpr contiguous_range<T>& operator= (std::array<U, N>& arr)
  {
    _start = arr.data();
    _size  = N;
    return *this;
  }

  template <
    class U,
    size_t N,
    std::enable_if_t<crange_type_is_const<U>>* = nullptr>
  constexpr contiguous_range<T>& operator= (std::array<U, N> const& arr)
  {
    _start = arr.data();
    _size  = N;
    return *this;
  }

  template <class U, size_t N>
  contiguous_range<T>& operator= (std::array<U, N>&& arr)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
    return *this;
  }

  // remember that contiguous_range is a non-owning reference that can get
  // easily invalidated...
  template <
    class U,
    class Alloc,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  contiguous_range (std::vector<U, Alloc>& vec)
    : contiguous_range {vec.data(), vec.size()}
  {}

  template <
    class U,
    class Alloc,
    std::enable_if_t<crange_type_is_const<U>>* = nullptr>
  contiguous_range (std::vector<U, Alloc> const& vec)
    : contiguous_range {vec.data(), vec.size()}
  {}

  template <class U, class Alloc>
  contiguous_range (std::vector<U, Alloc>&& vec)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
  }

  template <
    class U,
    class Alloc,
    std::enable_if_t<same_or_non_const_to_const<U>>* = nullptr>
  contiguous_range<T>& operator= (std::vector<U, Alloc>& vec)
  {
    _start = vec.data();
    _size  = vec.size();
    return *this;
  }

  template <
    class U,
    class Alloc,
    std::enable_if_t<crange_type_is_const<U>>* = nullptr>
  contiguous_range<T>& operator= (std::vector<U, Alloc>& vec)
  {
    _start = vec.data();
    _size  = vec.size();
    return *this;
  }

  template <class U, class Alloc>
  contiguous_range<T>& operator= (std::vector<U, Alloc>&& vec)
  {
    static_assert (!std::is_same_v<U, U>, "No binding to rvalues");
    return *this;
  }

  constexpr value_type& operator[] (size_t idx)
  {
    assert (_start);
    assert (idx < size());
    return _start[idx];
  }

  constexpr value_type const& operator[] (size_t idx) const
  {
    assert (_start);
    assert (idx < size());
    return _start[idx];
  }

  constexpr explicit operator bool() const { return !empty(); }

  // returns a copy/subrange with "count" elements dropped from the tail.
  constexpr contiguous_range<T> reduced (uint count) const
  {
    assert (count <= size());
    contiguous_range<T> r {*this};
    r._size -= count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }

  // returns a copy/subrange with "count" elements dropped from the head.
  constexpr contiguous_range<T> advanced (uint count) const
  {
    assert (count <= size());
    contiguous_range<T> r {*this};
    r._start += count;
    r._size -= count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  // returns a copy/subrange containing "count" elements starting from the
  // head.
  constexpr contiguous_range<T> get_head (uint count) const
  {
    assert (count <= size());
    contiguous_range<T> r {*this};
    r._size  = count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  // returns a copy/subrange containing "count" elements starting from the tail.
  constexpr contiguous_range<T> get_tail (uint count) const
  {
    assert (count <= size());
    contiguous_range<T> r {*this};
    r._start += r._size - count;
    r._size  = count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  // drops "count" elems from the head and returns the cut subrange
  constexpr contiguous_range<T> cut_head (uint count)
  {
    auto r = get_head (count);
    *this  = advanced (count);
    return r;
  }
  // drops "count" elems from the tail and returns the cut subrange
  constexpr contiguous_range<T> cut_tail (uint count)
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

  constexpr void clear()
  {
    _start = nullptr;
    _size == 0;
  }

  template <class U>
  contiguous_range<U> cast() const
  {
    constexpr auto big   = std::max (sizeof (value_type), sizeof (U));
    constexpr auto small = std::min (sizeof (value_type), sizeof (U));

    static_assert ((big % small) == 0, "sizes are not multiples");

    return {reinterpret_cast<U*> (_start), byte_size() / sizeof (U)};
  }
  // dummy parameter version to avoid calling make_crange(x).template cast<T>();
  template <class U>
  contiguous_range<U> cast (U) const
  {
    return cast<U>();
  }

  contiguous_range<std::add_const_t<T>> to_const() const { return *this; }

private:
  //----------------------------------------------------------------------------
  T*     _start = nullptr;
  size_t _size  = 0;
};
//------------------------------------------------------------------------------
template <class T>
static contiguous_range<T> make_contiguous_range (T* mem, size_t count)
{
  return {mem, count};
};

template <class T>
static contiguous_range<const T> make_contiguous_range (
  T const* mem,
  size_t   count)
{
  return {mem, count};
};

template <class T>
static contiguous_range<T> make_contiguous_range (T& mem)
{
  return {&mem, 1};
};

template <class T>
static contiguous_range<const T> make_contiguous_range (T const& mem)
{
  return {&mem, 1};
};

template <class T>
static void make_contiguous_range (T&& mem)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T, size_t N>
static contiguous_range<T> make_contiguous_range (
  T (&arr)[N],
  size_t count      = N,
  size_t offset_idx = 0)
{
  assert (count + offset_idx <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static contiguous_range<const T> make_contiguous_range (
  T const (&arr)[N],
  size_t count      = N,
  size_t offset_idx = 0)
{
  assert (count + offset_idx <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static contiguous_range<T> make_contiguous_range (
  std::array<T, N>& arr,
  size_t            count      = N,
  size_t            offset_idx = 0)
{
  assert ((count + offset_idx) <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static contiguous_range<const T> make_contiguous_range (
  std::array<T, N> const& arr,
  size_t                  count      = N,
  size_t                  offset_idx = 0)
{
  assert ((count + offset_idx) <= N && "out of bounds");
  return {arr.data() + offset_idx, count};
};

template <class T, size_t N>
static void make_contiguous_range (
  std::array<T, N>&& arr,
  size_t             count      = N,
  size_t             offset_idx = 0)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T, class Alloc>
static contiguous_range<T> make_contiguous_range (std::vector<T, Alloc>& vec)
{
  return {vec.data(), vec.size()};
};

template <class T, class Alloc>
static contiguous_range<const T> make_contiguous_range (
  std::vector<T, Alloc> const& vec)
{
  return {vec.data(), vec.size()};
};

template <class T, class Alloc>
static void make_contiguous_range (std::vector<T, Alloc>&& vec)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T, class Alloc>
static contiguous_range<T> make_contiguous_range (
  std::vector<T, Alloc>& vec,
  size_t                 count,
  size_t                 offset_idx = 0)
{
  assert (count + offset_idx <= vec.size() && "out of bounds");
  return {vec.data() + offset_idx, count};
};

template <class T, class Alloc>
static contiguous_range<const T> make_contiguous_range (
  std::vector<T, Alloc> const& vec,
  size_t                       count,
  size_t                       offset_idx = 0)
{
  assert (count + offset_idx <= vec.size() && "out of bounds");
  return {vec.data() + offset_idx, count};
};

template <class T, class Alloc>
static void make_contiguous_range (
  std::vector<T, Alloc>&& vec,
  size_t                  count,
  size_t                  offset_idx = 0)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T>
static contiguous_range<T> make_contiguous_range (
  contiguous_range<T> range,
  size_t              count,
  size_t              offset_idx = 0)
{
  assert (count + offset_idx <= range.size() && "out of bounds");
  return {&range[offset_idx], count};
};

template <class T>
static contiguous_range<const T> make_contiguous_range (
  contiguous_range<T> range,
  size_t              count,
  size_t              offset_idx = 0)
{
  assert (count + offset_idx <= range.size() && "out of bounds");
  return {&range[offset_idx], count};
};
//------------------------------------------------------------------------------
template <class T>
static void contiguous_range_memset (contiguous_range<T> range, int value = 0)
{
  memset (range.data(), value, range.size() * sizeof range[0]);
}
//------------------------------------------------------------------------------
template <class T, class U>
static uint contiguous_range_memcpy (
  contiguous_range<T>       dst,
  contiguous_range<const U> src)
{
  uint size = std::min (dst.size() * sizeof (T), src.size() * sizeof (U));
  memcpy (dst.data(), src.data(), size);
  return size;
}
// just for template deduction to work with const
template <class T, class U>
static uint contiguous_range_memcpy (
  contiguous_range<T> dst,
  contiguous_range<U> src)
{
  return contiguous_range_memcpy<T, const U> (dst, src);
}
//------------------------------------------------------------------------------
template <class T>
static uint contiguous_range_copy (
  contiguous_range<T>       dst,
  const contiguous_range<T> src)
{
  return contiguous_range_memcpy (dst, src) / sizeof (T);
}
//------------------------------------------------------------------------------
// shorter "contiguous_range" namings
//------------------------------------------------------------------------------
template <class T>
using crange = contiguous_range<T>;

template <class... Ts>
static auto make_crange (Ts&&... args)
{
  return make_contiguous_range (std::forward<Ts> (args)...);
}

template <class T>
static void crange_memset (crange<T> range, int value = 0)
{
  contiguous_range_memset (range, value);
}

template <class T, class U>
static uint crange_memcpy (crange<T> dst, const crange<U> src)
{
  return contiguous_range_memcpy (dst, src);
}

template <class T>
static uint crange_copy (crange<T> dst, const crange<T> src)
{
  return contiguous_range_copy (dst, src);
}
//------------------------------------------------------------------------------
} // namespace artv
