#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <type_traits>
#include <utility>
#include <vector>

#include <gcem.hpp>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

// TODO: Split this.

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

  constexpr contiguous_range<T> shrink_head (uint count) const
  {
    assert (count <= size());
    contiguous_range<T> r {*this};
    r._start += count;
    r._size -= count;
    return r;
  }

  constexpr contiguous_range<T> shrink_tail (uint count) const
  {
    assert (count <= size());
    contiguous_range<T> r {*this};
    r._size -= count;
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
template <class T>
static void contiguous_range_memset (crange<T> range, int value = 0)
{
  memset (range.data(), value, range.size() * sizeof range[0]);
}

template <class T>
static void crange_memset (crange<T> range, int value = 0)
{
  contiguous_range_memset (range, value);
}
//------------------------------------------------------------------------------
template <class T, class U>
static uint contiguous_range_memcpy (crange<T> dst, crange<const U> src)
{
  uint size = std::min (dst.size() * sizeof (T), src.size() * sizeof (U));
  memcpy (dst.data(), src.data(), size);
  return size;
}

template <class T, class U>
static uint crange_memcpy (crange<T> dst, crange<const U> src)
{
  return contiguous_range_memcpy (dst, src);
}
//------------------------------------------------------------------------------
template <class T>
static uint contiguous_range_copy (crange<T> dst, crange<const T> src)
{
  return crange_memcpy (dst, src) / sizeof (T);
}

template <class T>
static uint crange_copy (crange<T> dst, crange<const T> src)
{
  return crange_memcpy (dst, src) / sizeof (T);
}
//------------------------------------------------------------------------------
// Double pointers don't const convert, this is annoying on crange.
// https://stackoverflow.com/questions/5055655/double-pointer-const-correctness-warnings-in-c
template <class U, class T, size_t N>
static constexpr auto array_static_cast (std::array<T, N>& in)
{
  std::array<U, N> ret;
  for (uint i = 0; i < N; ++i) {
    ret[i] = static_cast<U> (in[i]);
  }
  return ret;
}

template <class U, class T, size_t N>
static constexpr auto array_const_cast (std::array<T, N>& in)
{
  std::array<U, N> ret;
  for (uint i = 0; i < N; ++i) {
    ret[i] = const_cast<U> (in[i]);
  }
  return ret;
}

template <class U, class T, size_t N>
static constexpr auto array_reinterpret_cast (std::array<T, N>& in)
{
  std::array<U, N> ret;
  for (uint i = 0; i < N; ++i) {
    ret[i] = reinterpret_cast<U> (in[i]);
  }
  return ret;
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
  constexpr T ln10 = gcem::log (10.);
  return db > m_inf_db ? std::exp (db * T {0.05} * ln10) : T {0.};
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T constexpr_gain_to_db (T gain, T m_inf_db = T {-130.})
{
  constexpr auto inv_ln10 = 1. / gcem::log (10.);
  return gain > T {0.} ? gcem::log (gain) * *T {20. * inv_ln10} : m_inf_db;
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
//------------------------------------------------------------------------------
// 0 to 44100/48000. 1 to 88/96KHz, etc... Assumes multiples of 44100 or 48000.
static uint get_samplerate_order (uint sample_rate)
{
  uint sr_order = sample_rate;
  if ((sr_order % 44100) != 0) {
    // assuming multiple of 48Khz
    auto srate_f = (double) sr_order;
    srate_f *= 44100. / 48000.;
    sr_order = (uint) srate_f;
    assert (sr_order % 44100 == 0 && "precission issues");
  }
  sr_order /= 44100; // 1, 2, 4, 8, 16, 32 ...
  sr_order = last_bit_set (sr_order); // 0, 1, 2, 3, 4, 5 ...
  return sr_order;
}
//------------------------------------------------------------------------------
template <uint N>
using k_uint = std::integral_constant<uint, N>;

template <int N>
using k_int = std::integral_constant<int, N>;

template <bool N>
using k_bool = std::integral_constant<bool, N>;

template <s8 N>
using k_s8 = std::integral_constant<s8, N>;

template <u8 N>
using k_u8 = std::integral_constant<u8, N>;

template <s16 N>
using k_s16 = std::integral_constant<s16, N>;

template <u16 N>
using k_u16 = std::integral_constant<u16, N>;

template <s32 N>
using k_s32 = std::integral_constant<s32, N>;

template <u32 N>
using k_u32 = std::integral_constant<u32, N>;

template <s64 N>
using k_s64 = std::integral_constant<s64, N>;

template <u64 N>
using k_u64 = std::integral_constant<u64, N>;
//------------------------------------------------------------------------------
namespace detail {
template <class tuple_like, class Func, size_t... Idxs>
auto tuple_unpack (
  tuple_like&& t,
  Func&&       unpack_f,
  std::index_sequence<Idxs...>)
{
  using ret_t = decltype (unpack_f (std::get<Idxs> (t)...));

  if constexpr (std::is_same_v<ret_t, void>) {
    unpack_f (std::get<Idxs> (t)...);
    return nullptr;
  }
  else {
    return unpack_f (std::get<Idxs> (t)...);
  }
}
} // namespace detail

template <template <class...> class tuple_like, class Func, class... Ts>
auto tuple_unpack (tuple_like<Ts...>&& t, Func&& unpack_f)
{
  return detail::tuple_unpack (
    std::forward<tuple_like<Ts...>> (t),
    std::forward<Func> (unpack_f),
    std::index_sequence_for<Ts...> {});
}
//------------------------------------------------------------------------------
template <class T>
constexpr bool is_aligned_to (uint align, T* v)
{
  assert (is_pow2 (align));
  auto addr = reinterpret_cast<same_size_uint<T*>> (v);
  return (addr & (align - 1)) == 0;
}

//------------------------------------------------------------------------------

} // namespace artv
