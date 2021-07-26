#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <type_traits>
#include <vector>

#include <gcem.hpp>

#include "artv-common/misc/short_ints.hpp"

namespace artv {

#ifdef __GNUC__
#define likely(x) __builtin_expect ((x), 1)
#define unlikely(x) __builtin_expect ((x), 0)
#else
#define likely(x) x
#define unlikely(x) x
#endif

#define VERSION_GET(major, minor, rev) (major * 1000000 + minor * 1000 + rev)

#define PP_EXPAND(x) x
#define PP_TO_STR(x) #x

// usage error_with_template_type<whatever> v; The compiler will print the type.
template <class T>
struct error_with_template_type;
//------------------------------------------------------------------------------
#define array_elems(x) sizeof (x) / sizeof (x[0])
//------------------------------------------------------------------------------
// mostly to be used on decltype statements.
template <class T, class... args>
static constexpr auto parampack_get_first (T v, args&&... vargs)
{
  return v;
}
//------------------------------------------------------------------------------
template <template <class...> class Dst, class Src>
struct pass_template_parameters;

template <
  template <class...>
  class Dst,
  template <class...>
  class Src,
  class... Ts>
struct pass_template_parameters<Dst, Src<Ts...>> {
  using type = Dst<Ts...>;
};
//------------------------------------------------------------------------------
template <class T>
static constexpr T* declptr()
{
  return (T*) 0;
}
//------------------------------------------------------------------------------
template <class T>
struct type_wrapper {
  using type = T;
};
//------------------------------------------------------------------------------
// unfortunately doesn't work with non type/class template parameters.
template <template <class...> class T, class U>
struct is_same_template : public std::false_type {};

template <template <class...> class T, class... Ts>
struct is_same_template<T, T<Ts...>> : public std::true_type {};

template <template <class...> class T, class... Ts>
struct is_same_template<T, const T<Ts...>> : public std::true_type {};

template <template <class...> class T, class U>
using is_same_template_t = typename is_same_template<T, U>::type;

template <template <class...> class T, class U>
static constexpr bool is_same_template_v = is_same_template<T, U>::value;

//------------------------------------------------------------------------------
template <class T = void, class... Ts>
static constexpr auto make_array (Ts&&... args)
{
  using U
    = std::conditional_t<std::is_same_v<T, void>, std::common_type_t<Ts...>, T>;
  static_assert (
    !std::disjunction_v<std::is_same<
      std::reference_wrapper<U>,
      std::remove_volatile_t<std::remove_const_t<Ts>>>...>,
    "reference_wrappers are not allowed");
  return std::array<U, sizeof...(Ts)> {std::forward<Ts> (args)...};
}

template <class... Ts>
static constexpr auto make_cstr_array (Ts&&... args)
{
  return make_array<char const*> (std::forward<Ts> (args)...);
}
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

  constexpr contiguous_range<T> shrink_head (uint count)
  {
    assert (count <= size());
    _start += count;
    _size -= count;
    return *this;
  }

  constexpr contiguous_range<T> shrink_tail (uint count)
  {
    assert (count <= size());
    _size -= count;
    return *this;
  }

  constexpr iterator       begin() { return _start; }
  constexpr const_iterator cbegin() { return _start; }
  constexpr const_iterator begin() const { return _start; }

  constexpr iterator       end() { return _start + _size; }
  constexpr const_iterator cend() { return _start + _size; }
  constexpr const_iterator end() const { return _start + _size; }

  constexpr size_t   size() const { return _size; }
  constexpr T*       data() { return _start; }
  constexpr T const* data() const { return _start; }

  constexpr bool empty() const { return size() == 0; }

private:
  //----------------------------------------------------------------------------
  T*     _start = nullptr;
  size_t _size  = 0;
};

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

template <class T>
static contiguous_range<T> make_contiguous_range (std::vector<T>& vec)
{
  return {vec.data(), vec.size()};
};

template <class T>
static contiguous_range<const T> make_contiguous_range (
  std::vector<T> const& vec)
{
  return {vec.data(), vec.size()};
};

template <class T>
static void make_contiguous_range (std::vector<T>&& vec)
{
  static_assert (!std::is_same_v<T, T>, "No binding to rvalues");
};

template <class T>
static contiguous_range<T> make_contiguous_range (
  std::vector<T>& vec,
  size_t          count,
  size_t          offset_idx = 0)
{
  assert (count + offset_idx <= vec.size() && "out of bounds");
  return {vec.data() + offset_idx, count};
};

template <class T>
static contiguous_range<const T> make_contiguous_range (
  std::vector<T> const& vec,
  size_t                count,
  size_t                offset_idx = 0)
{
  assert (count + offset_idx <= vec.size() && "out of bounds");
  return {vec.data() + offset_idx, count};
};

template <class T>
static void make_contiguous_range (
  std::vector<T>&& vec,
  size_t           count,
  size_t           offset_idx = 0)
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

// shorter "contiguous_range" naming -------------------------------------------
template <class T>
using crange = contiguous_range<T>;

template <class... Ts>
static auto make_crange (Ts&&... args)
{
  return make_contiguous_range (std::forward<Ts> (args)...);
}
//------------------------------------------------------------------------------
// Apply a function to each parameter on a  variadic argument pack
template <class funct>
static void apply (funct const& f)
{}

template <class funct, typename T, typename... args>
static void apply (funct const& f, T& v, args&&... vargs)
{
  f (v);
  apply (f, std::forward<args> (vargs)...);
}

template <class funct, typename T, typename... args>
static void apply (funct const& f, contiguous_range<T> r, args&&... vargs)
{
  for (int i = 0; i < r.size(); ++i) {
    f (r[i]);
  }
  apply (f, std::forward<args> (vargs)...);
}
//------------------------------------------------------------------------------
// divide memory on blocks. mainly to be used to wrap SIMD loop
template <class T, size_t N, class FunctorB, class FunctorT>
static void block_divide (
  uint              blocksize, // bytes
  std::array<T*, N> elems,
  uint              elem_count,
  FunctorB          f_blocks,
  FunctorT          f_tail)
{
  for (auto v : elems) {
    assert ((((uintptr_t) v) & (uintptr_t {blocksize} - 1)) == 0);
  }

  uint blocks = elem_count / (blocksize / sizeof (T));
  if (blocks) {
    f_blocks (elems, blocks);
  }
  uint offset = blocks * (blocksize / sizeof (T));
  uint tail   = elem_count - offset;
  if (tail) {
    for (auto& v : elems) {
      v += offset;
    }
    f_tail (elems, tail);
  }
}
//------------------------------------------------------------------------------
template <class T, size_t N, size_t Align>
class simd_mem {
public:
  using value_type                  = T;
  static constexpr size_t n_elems   = N;
  static constexpr size_t alignment = Align;

  T& operator[] (uint i)
  {
    assert (i < N);
    return mem[i];
  }
  T const& operator[] (uint i) const
  {
    assert (i < N);
    return mem[i];
  }
  T*               data() { return mem.data(); }
  T const*         data() const { return mem.data(); }
  constexpr size_t size() { return mem.size(); }

private:
  alignas (Align) std::array<T, N> mem;
};
//------------------------------------------------------------------------------
template <class T>
static constexpr T div_ceil (T num, T div)
{
  static_assert (std::is_integral<T>::value, "");
  return (num + div - 1) / div;
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T round_ceil (T num, T round)
{
  static_assert (std::is_integral<T>::value, "");
  return div_ceil (num, round) * round;
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T constexpr_db_to_gain (T db, T m_inf_db = T {-130.})
{
  return db > m_inf_db ? gcem::pow (T {10.0}, db * T {0.05}) : T {0.};
}
//------------------------------------------------------------------------------
template <class T>
static T db_to_gain (T db, T m_inf_db = T {-130.})
{
  return db > m_inf_db ? std::pow (T {10.0}, db * T {0.05}) : T {0.};
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T constexpr_gain_to_db (T gain, T m_inf_db = T {-130.})
{
  constexpr auto inv_log10 = 1. / gcem::log (10.);
  return gain > T {0.} ? gcem::log (gain) * *T {20. * inv_log10} : m_inf_db;
}
//------------------------------------------------------------------------------
template <class T>
static T gain_to_db (T gain, T m_inf_db = T {-130.})
{
  return gain > T {0.} ? std::log10 (gain) * T {20.} : m_inf_db;
}
//------------------------------------------------------------------------------
template <class T>
static T sgn_no_zero (T v, T neg = (T) -1., T pos_zero = (T) 1.)
{
  return (T) ((v < (T) 0) ? neg : pos_zero);
}
//------------------------------------------------------------------------------
template <class T>
static T sgn (T v)
{
  return (T) ((v > (T) 0) - (v < (T) 0));
}
//----------------------------------------------------------------------------
} // namespace artv
