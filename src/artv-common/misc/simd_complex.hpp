#pragma once

#include <array>
#include <cassert>
#include <complex>
#include <cstring>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <class V>
struct vec_complex {
  static_assert (is_vec_of_float_type_v<V>, "V must be a vector type.");
  //----------------------------------------------------------------------------
  static constexpr auto vec_size = vec_traits_t<V>::size;
  static constexpr auto size     = vec_size * 2; // total number of builtins
  using value_type               = V;
  using float_type               = vec_value_type_t<V>;
  using values_type = std::array<std::complex<float_type>, vec_size>;
  static_assert (sizeof (values_type) == (size * sizeof (float_type)), "?");
  //----------------------------------------------------------------------------
  V re;
  V im;
  //----------------------------------------------------------------------------
  constexpr vec_complex (V real = V {}, V imag = V {})
  {
    re = real;
    im = imag;
  }
  //----------------------------------------------------------------------------
  constexpr vec_complex (float_type real, float_type imag = float_type {})
  {
    re = vec_set<V> (real);
    im = vec_set<V> (imag);
  }
  //----------------------------------------------------------------------------
  constexpr vec_complex (vec_complex const& other)
  {
    re = other.re;
    im = other.im;
  }
  //----------------------------------------------------------------------------
  template <class U>
  constexpr vec_complex (vec_complex<U> const& other)
  {
    re = other.re;
    im = other.im;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator= (vec_complex const& other)
  {
    re = other.re;
    im = other.im;
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class U>
  vec_complex& operator= (vec_complex<U> const& other)
  {
    re = other.re;
    im = other.im;
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class U>
  vec_complex& operator= (value_type v)
  {
    re = v;
    im = vec_set<value_type> (float_type {});
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class U>
  vec_complex& operator= (float_type v)
  {
    re = vec_set<value_type> (v);
    im = vec_set<value_type> (float_type {});
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class U>
  vec_complex& operator+= (vec_complex<U> const& other)
  {
    re += other.re;
    im += other.im;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator+= (value_type v)
  {
    re += v;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator+= (float_type v)
  {
    *this += vec_set<value_type> (v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class U>
  vec_complex& operator-= (vec_complex<U> const& other)
  {
    re -= other.re;
    im -= other.im;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator-= (value_type v)
  {
    re -= v;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator-= (float_type v)
  {
    *this -= vec_set<value_type> (v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class U>
  vec_complex& operator*= (vec_complex<U> const& other)
  {
    // This is a naive implementation. The only reason it has an own
    // implementation is that the glibc version is equivalent, so at shuffling.
    // can be avoided<
    value_type re_new = re * other.re - im * other.im;
    im                = re * other.im + im * other.re;
    re                = re_new;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator*= (value_type v)
  {
    re *= v;
    im *= v;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator*= (float_type v)
  {
    *this *= vec_set<value_type> (v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class U>
  vec_complex& operator/= (vec_complex<U> const& other)
  {
    // This is a naive implementation. The only reason it has an own
    // implementation is that the glibc version is equivalent, so at shuffling.
    // can be avoided<
    value_type re_new = re * other.re + im * other.im;
    value_type den    = other.re * other.re + other.im * other.im;
    im                = (im * other.re - re * other.im) / den;
    re                = re_new / den;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator/= (value_type v)
  {
    re /= v;
    im /= v;
    return *this;
  }
  //----------------------------------------------------------------------------
  vec_complex& operator/= (float_type v)
  {
    *this /= vec_set<value_type> (v);
    return *this;
  }
  //----------------------------------------------------------------------------
  std::complex<float_type> to_std (uint idx) const
  {
    assert (idx < vec_size);
    return {re[idx], im[idx]};
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

template <class T>
struct is_complex_vec {
  static constexpr bool value = false;
};

template <class T>
struct is_complex_vec<vec_complex<T>> {
  static constexpr bool value = is_vec_v<T>;
};

template <class C>
static constexpr bool is_complex_vec_v = is_complex_vec<C>::value;

template <class C>
using enable_if_complex_vec_t = std::enable_if_t<is_complex_vec_v<C>>;

template <class C>
using complex_value_type_t = typename C::value_type;

template <class C>
using complex_float_type_t = typename C::float_type;
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> operator+ (vec_complex<V> const& v)
{
  return v;
}
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> operator- (vec_complex<V> const& v)
{
  return vec_complex<V> {-v.re, -v.im};
}
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> operator+ (
  vec_complex<V> const& lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r += rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator+ (vec_complex<V> const& lhs, V rhs)
{
  vec_complex<V> r {lhs};
  r += rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator+ (V lhs, vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r += rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator+ (
  vec_complex<V> const& lhs,
  vec_value_type_t<V>   rhs)
{
  vec_complex<V> r {lhs};
  r += rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator+ (
  vec_value_type_t<V>   lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r += rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> operator- (
  vec_complex<V> const& lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r -= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator- (vec_complex<V> const& lhs, V rhs)
{
  vec_complex<V> r {lhs};
  r -= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator- (V lhs, vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r -= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator- (
  vec_complex<V> const& lhs,
  vec_value_type_t<V>   rhs)
{
  vec_complex<V> r {lhs};
  r -= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator- (
  vec_value_type_t<V>   lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r -= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> operator* (
  vec_complex<V> const& lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r *= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator* (vec_complex<V> const& lhs, V rhs)
{
  vec_complex<V> r {lhs};
  r *= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator* (V lhs, vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r *= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator* (
  vec_complex<V> const& lhs,
  vec_value_type_t<V>   rhs)
{
  vec_complex<V> r {lhs};
  r *= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator* (
  vec_value_type_t<V>   lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r *= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> operator/ (
  vec_complex<V> const& lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r /= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator/ (vec_complex<V> const& lhs, V rhs)
{
  vec_complex<V> r {lhs};
  r /= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator/ (V lhs, vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r /= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator/ (
  vec_complex<V> const& lhs,
  vec_value_type_t<V>   rhs)
{
  vec_complex<V> r {lhs};
  r /= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> operator/ (
  vec_value_type_t<V>   lhs,
  vec_complex<V> const& rhs)
{
  vec_complex<V> r {lhs};
  r /= rhs;
  return r;
}
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> vec_polar (V radius, V theta = V {})
{
  return vec_complex<V> {radius * vec_cos (theta), radius * vec_sin (theta)};
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> vec_polar (
  vec_value_type_t<V> radius,
  vec_value_type_t<V> theta = vec_value_type_t<V> {})
{
  V re = vec_set<V> (radius * cos (theta));
  V im = vec_set<V> (radius * sin (theta));
  return vec_complex<V> {re, im};
}
//------------------------------------------------------------------------------
template <class V>
static vec_complex<V> vec_conj (vec_complex<V> v)
{
  return vec_complex<V> {v.re, -v.im};
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static vec_complex<V> vec_conj_mul (vec_complex<V> v)
{
  using T               = vec_value_type_t<V>;
  constexpr auto traits = vec_traits<V>();

  V re = vec_real (v);
  V im = vec_imag (v);
  re   = (re * re) + (im * im) - ((T) 2 * re) + (T) 1;
  im   = vec_set<V> ((T) 0);
  return vec_complex<V> {re, im};
}
//------------------------------------------------------------------------------
template <class V>
static V vec_real (vec_complex<V> v)
{
  return v.re;
}
//------------------------------------------------------------------------------
template <class V>
static V vec_imag (vec_complex<V> v)
{
  return v.im;
}
//------------------------------------------------------------------------------
// same interface as loading and storing from vectors. Basically wrapping memcpy
// and assertions.
//------------------------------------------------------------------------------
template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load (complex_float_type_t<C> const* src)
{
  C ret;
  memcpy (&ret, src, sizeof ret);
  return ret;
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load (C& dst, complex_float_type_t<C> const* src)
{
  dst = vec_load<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load (crange<const complex_float_type_t<C>> src)
{
  assert (src.size() >= C::size);
  return vec_load<C> (src.data());
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load (C& dst, crange<const complex_float_type_t<C>> src)
{
  dst = vec_load<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load (complex_value_type_t<C> const* src)
{
  C ret;
  memcpy (&ret, src, sizeof ret);
  return ret;
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load (C& dst, complex_value_type_t<C> const* src)
{
  dst = vec_load<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load (crange<const complex_value_type_t<C>> src)
{
  assert (src.size() >= C::size);
  return vec_load<C> (src.data());
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load (C& dst, crange<const complex_value_type_t<C>> src)
{
  dst = vec_load<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load_unaligned (complex_float_type_t<C> const* src)
{
  // dummy
  return vec_load<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load_unaligned (
  C&                             dst,
  complex_float_type_t<C> const* src)
{
  dst = vec_load_unaligned<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load_unaligned (crange<const complex_float_type_t<C>> src)
{
  assert (src.size() >= C::size);
  return vec_load_unaligned<C> (src.data());
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load_unaligned (
  C&                                    dst,
  crange<const complex_float_type_t<C>> src)
{
  dst = vec_load_unaligned<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load_unaligned (complex_value_type_t<C> const* src)
{
  // dummy
  return vec_load<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load_unaligned (
  C&                             dst,
  complex_value_type_t<C> const* src)
{
  dst = vec_load_unaligned<C> (src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline C vec_load_unaligned (crange<const complex_value_type_t<C>> src)
{
  assert (src.size() >= C::size);
  return vec_load_unaligned<C> (src.data());
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_load_unaligned (
  C&                                    dst,
  crange<const complex_value_type_t<C>> src)
{
  dst = vec_load_unaligned<C> (src);
}
//------------------------------------------------------------------------------
template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_store (complex_float_type_t<C>* dst, C src)
{
  memcpy (dst, &src, sizeof src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_store (crange<complex_float_type_t<C>> dst, C src)
{
  assert (dst.size() >= C::size);
  vec_store (dst.data(), src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_store (complex_value_type_t<C>* dst, C src)
{
  memcpy (dst, &src, sizeof src);
}

template <class C, enable_if_complex_vec_t<C>* = nullptr>
static inline void vec_store (crange<complex_value_type_t<C>> dst, C src)
{
  assert (dst.size() >= C::vec_size);
  vec_store (dst.data(), src);
}

//------------------------------------------------------------------------------
} // namespace artv
