// own old invention of (a dynamic extent) std::span, (previously named xspan)
// prepared from on C++20 probably inherit an std::span<T, std::dynamic_extent>
#pragma once

#include <array>
#include <iterator>
#include <type_traits>
#include <vector>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

template <class T>
class xspan;

namespace xspan_detail {

#if !defined(__cpp_lib_nonmember_container_access)
#error "This span implementation requires __cpp_lib_nonmember_container_access"
#endif

template <class T, class E>
static constexpr bool arrays_are_convertible_v = std::is_convertible_v<
  std::remove_pointer_t<decltype (std::data (std::declval<T>()))> (*)[],
  E (*)[]>;

template <class T>
static constexpr bool data_is_void_v = std::
  is_same_v<std::remove_cv_t<decltype (std::data (std::declval<T>()))>, void>;

template <class, class, class = void>
struct is_compatible : std::false_type {};

template <class T, class E>
struct is_compatible<
  T,
  E,
  std::enable_if_t<!data_is_void_v<T> && arrays_are_convertible_v<T, E>>>
  : std::true_type {};

template <class T, class E>
static constexpr bool is_compatible_v = is_compatible<T, E>::value;

template <class T, class E>
using enable_if_compatible = std::enable_if_t<is_compatible_v<T, E>>;
//------------------------------------------------------------------------------
template <class T>
struct is_std_array : public std::false_type {};

template <class T, size_t N>
struct is_std_array<std::array<T, N>> : public std::true_type {};

template <class T>
static constexpr bool is_std_array_v = is_std_array<T>::value;
//------------------------------------------------------------------------------
template <class T>
struct is_xspan : public std::false_type {};

template <class T>
struct is_xspan<xspan<T>> : public std::true_type {};

template <class T>
static constexpr bool is_xspan_v = is_xspan<T>::value;
//------------------------------------------------------------------------------
template <class, class = void>
struct callable_by_global_std_size_and_std_data : std::false_type {};

template <class T>
struct callable_by_global_std_size_and_std_data<
  T,
  std::void_t<
    decltype (std::size (std::declval<T>())),
    decltype (std::data (std::declval<T>()))>> : std::true_type {};

template <class T>
static constexpr bool callable_by_global_std_size_and_std_data_v
  = callable_by_global_std_size_and_std_data<T>::value;

template <typename T>
using uncvref_t = std::remove_cv_t<std::remove_reference<T>>;

// clang-format off
template <class T>
static constexpr bool is_container_v =
  !is_xspan_v<uncvref_t<T>> &&
  !std::is_array_v<uncvref_t<T>> &&
  !is_std_array_v<uncvref_t<T>> &&
  callable_by_global_std_size_and_std_data_v<T>;
// clang-format on

// For some reason in seems that the next template:
//
// template <class T, class E>
// using enable_if_compatible_container
//   = std::enable_if_t<is_container_v<T> && is_compatible_v<T, E>>;
//
// Is always evaluating both "is_container_v" and "is_compatible_v", when the
// expected behavior would be that if "is_container_v" is "false" then
// "is_compatible_v" shouldn't be evaluated.
//
// Without "container_compatible_workaround" every type fails on
// "data_is_void_v", as they have no "std::data" defined.
//
// This is on Clang14.
template <class T, class E, class = void>
struct container_compatible_workaround : std::false_type {};

template <class T, class E>
struct container_compatible_workaround<
  T,
  E,
  std::enable_if_t<is_container_v<T>>> : public is_compatible<T, E> {};

template <class T, class E>
using enable_if_compatible_container
  = std::enable_if_t<container_compatible_workaround<T, E>::value>;

//------------------------------------------------------------------------------
} // namespace xspan_detail
//------------------------------------------------------------------------------
// took this as reference:
// https://github.com/tcbrindle/span/blob/master/include/tcb/span.hpp
template <class T>
class xspan {
public:
  //----------------------------------------------------------------------------
  using element_type    = T;
  using value_type      = std::remove_cv_t<T>;
  using iterator        = element_type*;
  using const_iterator  = element_type const*;
  using pointer         = element_type*;
  using const_pointer   = element_type const*;
  using reference       = element_type&;
  using const_reference = element_type const&;
  using size_type       = size_t;

  static_assert (std::is_object_v<T>);
  static_assert (!std::is_abstract_v<T>);

  constexpr xspan()                         = default;
  ~xspan()                                  = default;
  constexpr xspan& operator= (xspan const&) = default;
  //----------------------------------------------------------------------------
  constexpr xspan (pointer ptr, size_type size) noexcept
  {
    _start = ptr;
    _size  = ptr ? size : 0;
  }
  //----------------------------------------------------------------------------
  constexpr xspan (pointer first, pointer last) noexcept
  {
    _start = first;
    assert (last >= first);
    _size = first - last;
  }
  //----------------------------------------------------------------------------
  template <
    size_type N,
    xspan_detail::
      enable_if_compatible<element_type (&)[N], element_type>* = nullptr>
  constexpr xspan (element_type (&arr)[N]) noexcept : xspan (arr, N)
  {}
  //----------------------------------------------------------------------------
  template <
    class U,
    size_type N,
    xspan_detail::
      enable_if_compatible<std::array<U, N>&, element_type>* = nullptr>
  constexpr xspan (std::array<U, N>& arr) noexcept : xspan (arr.data(), N)
  {}

  template <
    class U,
    size_type N,
    xspan_detail::
      enable_if_compatible<std::array<U, N> const&, element_type>* = nullptr>
  constexpr xspan (std::array<U, N> const& arr) noexcept : xspan (arr.data(), N)
  {}
  //----------------------------------------------------------------------------
  template <
    class Container,
    xspan_detail::
      enable_if_compatible_container<Container&, element_type>* = nullptr>
  constexpr xspan (Container& c) noexcept : xspan (std::data (c), std::size (c))
  {}

  template <
    class Container,
    xspan_detail::
      enable_if_compatible_container<Container const&, element_type>* = nullptr>
  constexpr xspan (Container const& c) noexcept
    : xspan (std::data (c), std::size (c))
  {}
  //----------------------------------------------------------------------------
  template <
    class U,
    std::enable_if_t<
      std::is_convertible_v<U (*)[], element_type (*)[]>>* = nullptr>
  constexpr xspan (xspan<U> const& s) noexcept : xspan (s.data(), s.size())
  {}
  //----------------------------------------------------------------------------
  constexpr reference at (size_type idx) noexcept
  {
    assert (_start);
    assert (idx < size());
    return _start[idx];
  }

  constexpr const_reference at (size_type idx) const noexcept
  {
    assert (_start);
    assert (idx < size());
    return _start[idx];
  }
  //----------------------------------------------------------------------------
  constexpr reference operator[] (size_type idx) noexcept { return at (idx); }
  constexpr const_reference operator[] (size_type idx) const noexcept
  {
    return at (idx);
  }
  //----------------------------------------------------------------------------
  constexpr explicit operator bool() const noexcept { return !empty(); }
  //----------------------------------------------------------------------------
  // returns a copy/subrange with "count" elements dropped from the tail.
  constexpr xspan<T> reduced (uint count) const noexcept
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._size -= count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  //----------------------------------------------------------------------------
  // returns a copy/subrange with "count" elements dropped from the head.
  constexpr xspan<T> advanced (uint count) const noexcept
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._start += count;
    r._size -= count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  //----------------------------------------------------------------------------
  // returns a copy/subrange containing "count" elements starting from the
  // head.
  constexpr xspan<T> get_head (uint count) const noexcept
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._size  = count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  //----------------------------------------------------------------------------
  // returns a copy/subrange containing "count" elements starting from the tail.
  constexpr xspan<T> get_tail (uint count) const noexcept
  {
    assert (count <= size());
    xspan<T> r {*this};
    r._start += r._size - count;
    r._size  = count;
    r._start = r._size ? r._start : nullptr;
    return r;
  }
  //----------------------------------------------------------------------------
  // drops "count" elems from the head and returns the cut subrange
  constexpr xspan<T> cut_head (uint count) noexcept
  {
    auto r = get_head (count);
    *this  = advanced (count);
    return r;
  }
  //----------------------------------------------------------------------------
  // drops "count" elems from the tail and returns the cut subrange
  constexpr xspan<T> cut_tail (uint count) noexcept
  {
    auto r = get_tail (count);
    *this  = reduced (count);
    return r;
  }
  //----------------------------------------------------------------------------
  constexpr iterator       begin() noexcept { return _start; }
  constexpr const_iterator cbegin() noexcept { return _start; }
  constexpr const_iterator begin() const noexcept { return _start; }
  //----------------------------------------------------------------------------
  constexpr iterator       end() noexcept { return _start + _size; }
  constexpr const_iterator cend() noexcept { return _start + _size; }
  constexpr const_iterator end() const noexcept { return _start + _size; }
  //----------------------------------------------------------------------------
  constexpr size_type size() const noexcept { return _size; }
  constexpr size_type byte_size() const noexcept
  {
    return _size * sizeof (element_type);
  }
  constexpr pointer       data() noexcept { return _start; }
  constexpr const_pointer data() const noexcept { return _start; }
  //----------------------------------------------------------------------------
  constexpr bool empty() const noexcept { return size() == 0; }
  //----------------------------------------------------------------------------
  constexpr reference       front() noexcept { return at (0); }
  constexpr const_reference front() const noexcept { return at (0); }
  constexpr reference       back() noexcept { return at (_size - 1); }
  constexpr const_reference back() const noexcept { return at (_size - 1); }
  //----------------------------------------------------------------------------
  constexpr auto first (size_type size) const noexcept
  {
    return get_head (size);
  }
  constexpr auto last (size_type size) const noexcept
  {
    return get_tail (size);
  }
  //----------------------------------------------------------------------------
  constexpr auto subspan (size_type offset) const { return advanced (offset); }
  constexpr auto subspan (size_type offset, size_type count) const
  {
    return advanced (offset).get_head (count);
  }
  //----------------------------------------------------------------------------
  constexpr void clear() noexcept
  {
    _start = nullptr;
    _size == 0;
  }
  //----------------------------------------------------------------------------
  template <class U>
  xspan<U> cast() const noexcept
  {
    constexpr auto big   = std::max (sizeof (element_type), sizeof (U));
    constexpr auto small = std::min (sizeof (element_type), sizeof (U));

    static_assert ((big % small) == 0, "sizes are not multiples");

    return {reinterpret_cast<U*> (_start), byte_size() / sizeof (U)};
  }
  //----------------------------------------------------------------------------
  // dummy parameter version to avoid calling xspan(x).template cast<T>();
  template <class U>
  xspan<U> cast (U) const noexcept
  {
    return cast<U>();
  }
  //----------------------------------------------------------------------------
  xspan<std::add_const_t<T>> to_const() const noexcept { return *this; }

private:
  //----------------------------------------------------------------------------
  T*        _start = nullptr;
  size_type _size  = 0;
};
//------------------------------------------------------------------------------
template <class T>
static void xspan_memset (xspan<T> range, int value = 0)
{
  memset (range.data(), value, range.size() * sizeof range[0]);
}
//------------------------------------------------------------------------------
template <class T, class U>
static uint xspan_memcpy (xspan<T> dst, xspan<U const> src)
{
  uint size = std::min (dst.size() * sizeof (T), src.size() * sizeof (U));
  memcpy (dst.data(), src.data(), size);
  return size;
}
// just for template deduction to work with const
template <class T, class U>
static uint xspan_memcpy (xspan<T> dst, xspan<U> src)
{
  return xspan_memcpy<T, U const> (dst, src);
}
//------------------------------------------------------------------------------
// dangerous variant, but not more than memcpy itself, which I'm already using.
template <class T, class U>
static uint xspan_memdump (T* dst, xspan<U const> src)
{
  uint size = src.size() * sizeof (U);
  memcpy (dst, src.data(), src.size() * sizeof (U));
  return size;
}
// just for template deduction to work with const
template <class T, class U>
static uint xspan_memdump (T* dst, xspan<U> src)
{
  return xspan_memdump<T, U const> (dst, src);
}
//------------------------------------------------------------------------------
template <class T>
static uint xspan_copy (xspan<T> dst, xspan<T> const src)
{
  return xspan_memcpy (dst, src) / sizeof (T);
}
//------------------------------------------------------------------------------
#if defined(__cpp_deduction_guides)

template <class T, size_t N>
xspan (T (&)[N]) -> xspan<T>;

template <class T, size_t N>
xspan (std::array<T, N>&) -> xspan<T>;

template <class T, size_t N>
xspan (std::array<T, N> const&) -> xspan<T const>;

template <class Container>
xspan (Container&) -> xspan<
  std::remove_reference_t<decltype (*std::data (std::declval<Container&>()))>>;

template <class Container>
xspan (Container const&) -> xspan<const typename Container::value_type>;

#else
#error "This span implementation requires __cpp_deduction_guides"
#endif
} // namespace artv
