#pragma once

#include <algorithm>
#include <assert.h>
#include <gcem.hpp>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

namespace detail {

template <uint N_int, uint N_frac, class = void>
struct fp_type {
  static_assert (
    (N_int + N_frac) <= 64,
    "only types representable by 64 bits supported");
};

template <uint N_int, uint N_frac>
struct fp_type<
  N_int,
  N_frac,
  std::enable_if_t<
    std::clamp<uint> (N_int + N_frac, 33, 64) == (N_int + N_frac)>> {
  using type = u64;
};

template <uint N_int, uint N_frac>
struct fp_type<
  N_int,
  N_frac,
  std::enable_if_t<
    std::clamp<uint> (N_int + N_frac, 17, 32) == (N_int + N_frac)>> {
  using type = u32;
};

template <uint N_int, uint N_frac>
struct fp_type<
  N_int,
  N_frac,
  std::enable_if_t<
    std::clamp<uint> (N_int + N_frac, 9, 16) == (N_int + N_frac)>> {
  using type = u16;
};

template <uint N_int, uint N_frac>
struct fp_type<
  N_int,
  N_frac,
  std::enable_if_t<
    std::clamp<uint> (N_int + N_frac, 0, 8) == (N_int + N_frac)>> {
  using type = u8;
};

} // namespace detail

//------------------------------------------------------------------------------
// A (not that) quick and dirty class to avoid some mistakes when doing fixed
// point arithmetic. It could probably be good to use a 3rd party lib as there
// are a lot of details and decisions to be made.
//
// The design criteria is that every constructor and operator overload including
// assignments allows making no mistakes on the ranges.
//------------------------------------------------------------------------------
template <uint N_int, uint N_frac, bool is_signed = true>
class fixed_point {
public:
  // incomplete, expand as necessary.
  static constexpr uint n_int  = N_int;
  static constexpr uint n_frac = N_frac;

  using value_type_uint = typename detail::fp_type<n_int, n_frac>::type;

  using value_type = std::conditional_t<
    is_signed,
    value_type_uint,
    std::make_signed_t<value_type_uint>>;

  using float_type
    = std::conditional_t<sizeof (value_type) <= 2, float, double>;

  static constexpr auto float_factor
    = (float_type) gcem::pow (2., (float_type) n_frac);

  static constexpr value_type int_max
    = lsb_mask<value_type_uint> (n_int - (!!n_int & is_signed));

  static constexpr value_type int_min = (is_signed && n_int)
    ? (value_type) -lsb_mask<value_type_uint> (n_int - 1)
    : 0;

  static constexpr float_type flt_max
    = ((float_type) int_max) + (float_type) 1 - (float_type) 1 / float_factor;

  static constexpr float_type flt_min
    = ((float_type) int_min) - (float_type) 1 + (float_type) 1 / float_factor;
  //----------------------------------------------------------------------------
  fixed_point() { _v = decltype (_v) {}; }
  //----------------------------------------------------------------------------
  template <
    uint N_int_v,
    uint N_frac_v,
    std::enable_if_t<N_int_v <= n_int && N_frac_v <= n_frac>* = nullptr>
  constexpr fixed_point (fixed_point<N_int_v, N_frac_v> fp)
  {
    *this = fp;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  constexpr fixed_point (T flt)
  {
    *this = flt;
  }
  //----------------------------------------------------------------------------
  template <
    class T,
    std::enable_if_t<
      std::is_integral_v<T> && std::is_signed_v<T> == is_signed>* = nullptr>
  constexpr fixed_point (T intv)
  {
    *this = intv;
  }
  //----------------------------------------------------------------------------
  ~fixed_point() = default;
  //----------------------------------------------------------------------------
  template <
    uint N_int_v,
    uint N_frac_v,
    std::enable_if_t<N_int_v <= n_int && N_frac_v <= n_frac>* = nullptr>
  constexpr fixed_point& operator= (fixed_point<N_int_v, N_frac_v> fp)
  {
    constexpr auto shift = (int) (N_frac_v - n_frac);

    if constexpr (shift > 0) {
      // source has more fractional bits
      _v = fp._v >> shift;
    }
    else {
      // destination has more or equal fractional bits
      _v = fp._v;
      _v <<= -shift;
    }
    return *this;
  }

  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  constexpr fixed_point& operator= (T flt)
  {
    // Assertions probably not 100% correct/verified until the corner case of
    // the last decimals betweem the fractional and integral part but enough to
    // catch mistakes.
    assert (flt <= flt_max);
    assert (flt >= flt_min);

    _v = (value_type) (flt * float_factor);

    return *this;
  }
  //----------------------------------------------------------------------------
  template <
    class T,
    std::enable_if_t<
      std::is_integral_v<T> && std::is_signed_v<T> == is_signed>* = nullptr>
  constexpr fixed_point& operator= (T intv)
  {
    assert (intv <= int_max);
    assert (intv >= int_min);

    if constexpr (sizeof intv > sizeof _v) {
      _v = (value_type) (intv << n_frac);
    }
    else {
      _v = ((value_type) intv) << n_frac;
    }
    return *this;
  }
  //----------------------------------------------------------------------------
  // All the operator overloads are lossless
  //----------------------------------------------------------------------------
  template <
    uint N_int_v,
    uint N_frac_v,
    std::enable_if_t<N_int_v <= n_int && N_frac_v <= n_frac>* = nullptr>
  constexpr auto operator+= (fixed_point<N_int_v, N_frac_v, is_signed> rhs)
  {
    _v += ((value_type) rhs._v) << (n_frac - N_frac_v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator+ (fixed_point<N_int_v, N_frac_v, is_signed> rhs)
  {
    constexpr uint op_n_int  = std::max (n_int, N_int_v);
    constexpr uint op_n_frac = std::max (n_frac, N_frac_v);

    using lossless_type = fixed_point<op_n_int, op_n_frac>;
    using value_type    = typename lossless_type::value_type;

    constexpr auto shift = (int) (N_frac_v - n_frac);

    lossless_type ret;
    if constexpr (shift > 0) {
      // rhs has more fractional bits
      ret._v = rhs._v;
      ret._v += ((value_type) _v) << shift;
    }
    else {
      // instance has more or equal fractional bits
      ret._v = _v;
      ret._v += ((value_type) rhs._v) << -shift;
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator* (fixed_point<N_int_v, N_frac_v, is_signed> rhs)
  {
    constexpr uint op_n_int  = n_int + N_int_v;
    constexpr uint op_n_frac = n_frac + N_frac_v;

    using lossless_type = fixed_point<op_n_int, op_n_frac>;
    using value_type    = typename lossless_type::value_type;

    lossless_type ret;
    ret._v = _v * rhs._v;
    return ret;
  }
  //----------------------------------------------------------------------------
  // The operators with no overload keep the type unchanged, they only allow
  // to operate with types of less or equal resolution on the integer part.
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr void mul (fixed_point<N_int_v, N_frac_v, is_signed> rhs)
  {
    static_assert (
      (n_int + n_frac + N_frac_v) <= (sizeof _v * 8),
      "can't keep the required width for the integer part");
    _v *= rhs._v;
    _v >>= N_frac_v;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr void add (fixed_point<N_int_v, N_frac_v, is_signed> rhs)
  {
    static_assert (n_int >= N_int_v, "Integer resolution of rhs is bigger");

    constexpr auto shift = (int) (N_frac_v - n_frac);

    if constexpr (shift > 0) {
      // rhs has more fractional bits
      _v += (value_type) (rhs._v >> shift);
    }
    else {
      // instance has more or equal fractional bits
      _v += ((value_type) rhs._v) << -shift;
    }
  }
  //----------------------------------------------------------------------------
  constexpr value_type integer() const
  {
    return (_v >> n_frac) & lsb_mask<value_type_uint> (n_int);
  }
  //----------------------------------------------------------------------------
  template <class T = std::conditional_t<n_frac <= 23, float, double>>
  constexpr auto fraction() const
  {
    constexpr double factor = (T) 1. / float_factor;
    return factor * (T) (_v & lsb_mask<value_type_uint> (n_frac));
  }
  //----------------------------------------------------------------------------
  constexpr void       set_raw (value_type v) { _v = v; }
  constexpr value_type get_raw() const { return _v; }
  //----------------------------------------------------------------------------
private:
  template <uint, uint, bool>
  friend class fixed_point;

  value_type _v;
};
//------------------------------------------------------------------------------
}; // namespace artv
