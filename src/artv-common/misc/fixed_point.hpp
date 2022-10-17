#pragma once

#include <algorithm>
#include <assert.h>
#include <gcem.hpp>
#include <type_traits>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <class S, class U, class F, uint N_bits = sizeof (S) * 8>
struct fp_int_promotion_step {
  static_assert (sizeof (S) == sizeof (U));
  using tsigned                  = S;
  using tunsigned                = U;
  using tfloat                   = F;
  static constexpr size_t n_bits = N_bits; // might be a vector type
};
//------------------------------------------------------------------------------
struct std_fp_types_trait {
  using promotions = mp_list<
    fp_int_promotion_step<s8, u8, float>,
    fp_int_promotion_step<s16, u16, float>,
    fp_int_promotion_step<s32, u32, double>,
    fp_int_promotion_step<s64, u64, long double>
    //  fp_int_promotion_step<__int128, unsigned __int128, long double>
    >;
  static constexpr uint vector_size = 0;
};
//------------------------------------------------------------------------------
template <class T_signed, class T_unsigned, class T_float, uint Vec_Size = 0>
struct single_fp_type_trait {
  using promotions
    = mp_list<fp_int_promotion_step<T_signed, T_unsigned, T_float>>;
  static constexpr uint vector_size = Vec_Size;
};
//------------------------------------------------------------------------------
namespace detail {
//------------------------------------------------------------------------------
template <uint N>
struct first_type_n_bits_gt_n {
  template <class U>
  using fn = mp11::mp_bool<U::n_bits >= N>;
};
//------------------------------------------------------------------------------
template <
  uint N_sign,
  uint N_int,
  uint N_frac,
  class Traits = std_fp_types_trait>
class fp_int_base {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_int     = N_int;
  static constexpr uint n_frac    = N_frac;
  static constexpr uint n_sign    = N_sign;
  static constexpr uint n_bits    = n_sign + n_int + n_frac;
  static constexpr uint is_signed = n_sign == 1;
  using traits                    = Traits;
  using promotions                = typename Traits::promotions;
  //----------------------------------------------------------------------------
  using type_idx
    = mp11::mp_find_if_q<promotions, first_type_n_bits_gt_n<n_bits>>;

  static_assert (n_sign <= 1, "Invalid value. 1 is signed, 0 unsigned");
  static_assert (
    type_idx::value < mp11::mp_size<promotions>::value,
    "Number of bits required not representable by any available type");
  //----------------------------------------------------------------------------
  using type          = mp11::mp_at<promotions, type_idx>;
  using unsigned_type = typename type::tunsigned;
  using signed_type   = typename type::tsigned;
  using float_type    = typename type::tfloat;
  using value_type = std::conditional_t<is_signed, signed_type, unsigned_type>;
  //----------------------------------------------------------------------------
  // can't use "std::numeric_limits" or assume constructors based on a float
  // type. These look that weird because they work with Clang and Gcc Types

  static constexpr unsigned_type unsigned_zeroes = unsigned_type {};
  static constexpr unsigned_type unsigned_ones   = unsigned_type {} - 1U;
  static constexpr unsigned_type unsigned_one    = unsigned_type {} + 1U;

  static constexpr auto float_factor = (float_type) (unsigned_one << n_frac);

  static constexpr value_type int_max
    = (value_type) (~(unsigned_ones << (n_int + n_sign)));

  static constexpr value_type int_min
    = is_signed ? -int_max - 1 : value_type {};

  // these limits might not always be precise depending on the floating point
  // type
  static constexpr float_type float_one  = float_type {} + 1.f;
  static constexpr float_type float_imin = float_type {} + int_min;
  static constexpr float_type float_imax = float_type {} + int_max;

  static constexpr float_type flt_max
    = float_imax + float_one - float_one / float_factor;

  static constexpr float_type flt_min = (is_signed)
    ? float_imin - float_one + float_one / float_factor
    : float_type {};

  //----------------------------------------------------------------------------
  template <class T>
  using enable_if_int_same_signedness = std::enable_if_t<
    std::is_integral_v<T> && std::is_signed_v<T> == (bool) n_sign>;
  //----------------------------------------------------------------------------
  fp_int_base() { _v = decltype (_v) {}; }
  //----------------------------------------------------------------------------
  ~fp_int_base() = default;
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr fp_int_base& operator= (
    fp_int_base<n_sign, N_int_v, N_frac_v, promotions> fp)
  {
    constexpr auto shift = (int) (N_frac_v - n_frac);

    if constexpr (shift > 0) {
      // source has more fractional bits
      _v = ashr<shift> (fp._v);
    }
    else {
      // destination has more or equal fractional bits
      _v = fp._v;
      _v = ashl<-shift> (fp._v);
    }
    return *this;
  }
  //----------------------------------------------------------------------------
  static constexpr fp_int_base from (value_type raw_val)
  {
    fp_int_base r;
    r._v = raw_val;
    return r;
  }
  //----------------------------------------------------------------------------
  // lossy conversion
  static constexpr fp_int_base from_float (float_type flt)
  {
    // Assertions probably not 100% correct/verified, as the non-constant
    // float/double resolution makes this a non-trivial test. For the meantime
    // this might be enough to catch most mistakes. TODO
    auto fx = flt_max;
    auto fn = flt_min;
    auto ff = float_factor;

    if constexpr (traits::vector_size == 0) {
      assert (flt <= flt_max);
      assert (flt >= flt_min);
    }
    else {
      auto lte_max = flt <= flt_max;
      auto gte_min = flt >= flt_min;
      for (uint i = 0; i < traits::vector_size; ++i) {
        assert (lte_max[i]);
        assert (gte_min[i]);
      }
    }

    fp_int_base r;
    float_type  v = flt * float_factor;
    v += float_type {
      (v > 0) ? float_type {0.5} : float_type {-0.5}}; // rounding
    r._v = (value_type) v;
    return r;
  }
  //----------------------------------------------------------------------------
  // value_type is taken as an integer with no scaling.
  static constexpr fp_int_base from_int (value_type intv)
  {
    assert (intv <= int_max);
    assert (intv >= int_min);

    if constexpr (traits::vector_size == 0) {
      assert (intv <= int_max);
      assert (intv >= int_min);
    }
    else {
      auto lte_max = intv <= int_max;
      auto gte_min = intv >= int_min;
      for (uint i = 0; i < traits::vector_size; ++i) {
        assert (lte_max[i]);
        assert (gte_min[i]);
      }
    }

    fp_int_base r;
    if constexpr (sizeof intv > sizeof _v) {
      r._v = (value_type) ashl<n_frac> (intv);
    }
    else {
      r._v = ashl<n_frac> ((value_type) intv);
    }
    return r;
  }
  //----------------------------------------------------------------------------
  constexpr void load (value_type v) { *this = from (v); }
  //----------------------------------------------------------------------------
  constexpr void load_float (float_type v) { *this = from_float (v); }
  //----------------------------------------------------------------------------
  constexpr void load_int (value_type v) { *this = from_int (v); }
  //----------------------------------------------------------------------------
  template <class T>
  constexpr T cast (T) const
  {
    constexpr auto shift = (int) (T::n_frac - n_frac);
    T              ret;
    using int_type = typename T::value_type;
    if constexpr (shift > 0) {
      // dst has more fractional bits
      ret._v = ashl<shift> ((int_type) _v);
    }
    else {
      // instance has more or equal fractional bits
      ret._v = (int_type) ashr<-shift> (_v);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T>
  constexpr T cast() const
  {
    return cast (T {});
  }
  //----------------------------------------------------------------------------
  constexpr float_type as_float() const
  {
    constexpr float_type factor = (float_type) 1. / float_factor;
    return factor * (float_type) (_v);
  }
  //----------------------------------------------------------------------------
  constexpr value_type as_int() const { return (_v >> n_frac); }
  //----------------------------------------------------------------------------
  constexpr float_type fraction() const
  {
    constexpr double factor = (float_type) 1. / float_factor;
    return factor * (float_type) (_v & lsb_mask<unsigned_type> (n_frac));
  }
  //----------------------------------------------------------------------------
  // useful E.g. after random loads to check that the value is in range. Ideally
  // should happen after each load, store, constructor, assignment, etc from
  // external sources, unfortunately this is a numeric/performace sensitive
  // class.
  void normalize() { _v = alsb_mask (_v, n_frac + n_int); }
  //----------------------------------------------------------------------------
  // TODO: comparison operators
  //----------------------------------------------------------------------------
  constexpr value_type value() const { return _v; }
  //----------------------------------------------------------------------------
  template <uint, uint, uint, class>
  friend class fp_int_base;

  value_type _v;
};
//------------------------------------------------------------------------------

} // namespace detail
//------------------------------------------------------------------------------
template <
  uint N_sign,
  uint N_int,
  uint N_frac,
  bool Is_lossless = true,
  class Traits     = std_fp_types_trait>
class fixed_point_int;

template <class T>
struct is_fixed_point_int : public std::false_type {};

template <uint N_sign, uint N_int, uint N_frac, bool Is_lossless, class Traits>
struct is_fixed_point_int<
  fixed_point_int<N_sign, N_int, N_frac, Is_lossless, Traits>>
  : public std::true_type {};

template <class T>
static constexpr bool is_fixed_point_int_v = is_fixed_point_int<T>::value;

//------------------------------------------------------------------------------
// A simple fixed point arithmetic class that always stays with the same
// precision and bit width (lossy).
//------------------------------------------------------------------------------
template <uint N_sign, uint N_int, uint N_frac, class Traits>
class fixed_point_int<N_sign, N_int, N_frac, false, Traits>
  : private detail::fp_int_base<N_sign, N_int, N_frac, Traits> {
private:
  using base = detail::fp_int_base<N_sign, N_int, N_frac, Traits>;

public:
  //----------------------------------------------------------------------------
  using float_type = typename base::float_type;
  using value_type = typename base::value_type;
  using traits     = typename base::traits;
  using promotions = typename base::promotions;
  //----------------------------------------------------------------------------
  static constexpr bool       is_lossless  = false;
  static constexpr auto       n_bits       = base::n_bits;
  static constexpr auto       n_int        = base::n_int;
  static constexpr auto       n_frac       = base::n_frac;
  static constexpr auto       n_sign       = base::n_sign;
  static constexpr auto       is_signed    = base::is_signed;
  static constexpr auto       float_factor = base::float_factor;
  static constexpr value_type int_max      = base::int_max;
  static constexpr value_type int_min      = base::int_min;
  static constexpr float_type flt_max      = base::flt_max;
  static constexpr float_type flt_min      = base::flt_min;
  //----------------------------------------------------------------------------
  using base::as_float;
  using base::as_int;
  using base::base;
  using base::fraction;
  using base::load;
  using base::load_float;
  using base::load_int;
  using base::normalize;
  using base::value;
  using base::operator=;
  //----------------------------------------------------------------------------
  template <uint N_intb, uint N_fracb>
  using compatiblefp = fixed_point_int<n_sign, N_intb, N_fracb, false, traits>;
  //----------------------------------------------------------------------------
  // raw constructor
  explicit fixed_point_int (value_type v) { this->_v = v; }
  //----------------------------------------------------------------------------
  template <
    uint N_int_v,
    uint N_frac_v,
    std::enable_if_t<N_int_v <= n_int && N_frac_v <= n_frac>* = nullptr>
  constexpr fixed_point_int (compatiblefp<N_int_v, N_frac_v> fp)
  {
    *this = fp;
  }
  //----------------------------------------------------------------------------
  ~fixed_point_int() = default;
  //----------------------------------------------------------------------------
  // All the operators might lose precission
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr fixed_point_int& operator= (compatiblefp<N_int_v, N_frac_v> fp)
  {
    *to_base (this) = *to_base (&fp);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto& operator+= (compatiblefp<N_int_v, N_frac_v> rhs)
  {
    constexpr auto shift = (int) (N_frac_v - n_frac);
    if constexpr (shift > 0) {
      this->_v += (value_type) ashr<shift> (rhs._v);
    }
    else {
      this->_v += ashl<-shift> ((value_type) rhs._v);
    }
    return *this;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator+ (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    auto ret {*this};
    ret += rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto& operator-= (compatiblefp<N_int_v, N_frac_v> rhs)
  {
    constexpr auto shift = (int) (N_frac_v - n_frac);
    if constexpr (shift > 0) {
      this->_v -= (value_type) ashr<shift> (rhs._v);
    }
    else {
      this->_v -= ashl<-shift> ((value_type) rhs._v);
    }
    return *this;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator- (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    auto ret {*this};
    ret -= rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto& operator*= (compatiblefp<N_int_v, N_frac_v> rhs)
  {
    // only make space to handle (and revert) the shift caused by multiplying.
    constexpr uint op_n_int  = n_int + N_frac_v;
    constexpr uint op_n_frac = n_frac;

    using dst_type = compatiblefp<op_n_int, op_n_frac>;
    using bigger   = typename dst_type::value_type;
    this->_v       = ashr<N_frac_v> ((bigger) this->_v * (bigger) rhs._v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator* (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    auto ret {*this};
    ret *= rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto& operator/= (compatiblefp<N_int_v, N_frac_v> rhs)
  {
    // only make space to handle the pre-shift required for dividing while
    // keeping the same scaling.
    constexpr uint op_n_int  = n_int + N_frac_v;
    constexpr uint op_n_frac = n_frac;
    using dst_type           = fixed_point_int<n_sign, op_n_int, op_n_frac>;

    auto tmp = cast<dst_type>();
    this->_v = (value_type) (ashl<N_frac_v> (tmp._v) / rhs._v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator/ (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    auto ret {*this};
    ret /= rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T>
  constexpr T cast (T) const
  {
    static_assert (is_fixed_point_int_v<T>);
    static_assert (T::is_signed == is_signed);
    static_assert (std::is_same_v<typename T::promotions, promotions>);

    using ptr       = T*;
    using base_type = std::remove_pointer_t<decltype (to_base (ptr {nullptr}))>;
    T ret;
    ret._v = base::template cast<base_type>()._v;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T>
  constexpr T cast() const
  {
    return cast (T {});
  }
  //----------------------------------------------------------------------------
  static constexpr fixed_point_int from (value_type raw_val)
  {
    return fixed_point_int {raw_val};
  }
  //----------------------------------------------------------------------------
  static constexpr fixed_point_int from_float (float_type flt)
  {
    return fixed_point_int {base::from_float (flt)._v};
  }
  //----------------------------------------------------------------------------
  static constexpr fixed_point_int from_int (value_type intv)
  {
    return fixed_point_int {base::from_int (intv)._v};
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <uint S, uint I, uint F>
  static constexpr auto to_base (fixed_point_int<S, I, F, false, traits>* v)
  {
    return static_cast<detail::fp_int_base<S, I, F>*> (v);
  }
  //----------------------------------------------------------------------------
  template <uint, uint, uint, bool, class>
  friend class fixed_point_int;
};
//------------------------------------------------------------------------------
// A simple fixed point arithmetic class that dynamically expands its range
// fixed point and range depending on the operations.
//
// Notice that it doesn't follow the ranges over what the maximum that a builtin
// type can represent.
//
// With other wording, this doesn't allow the scaling (fixed point) to be over
// 2^(biggest type bits), while it conceptualy would be perfectly possible with
// another kind of abstraction, e.g. passing a precision, a builtin type and a
// factor (fixed point location).
//
// Anyhow, for that it would probably be better to bite the bullet and add a
// dependency on some mature existing library.
//------------------------------------------------------------------------------
template <uint N_sign, uint N_int, uint N_frac, class Traits>
class fixed_point_int<N_sign, N_int, N_frac, true, Traits>
  : private detail::fp_int_base<N_sign, N_int, N_frac, Traits> {
private:
  using base = detail::fp_int_base<N_sign, N_int, N_frac, Traits>;

public:
  //----------------------------------------------------------------------------
  using float_type = typename base::float_type;
  using value_type = typename base::value_type;
  using traits     = typename base::traits;
  using promotions = typename base::promotions;
  //----------------------------------------------------------------------------
  static constexpr bool       is_lossless  = true;
  static constexpr auto       n_bits       = base::n_bits;
  static constexpr auto       n_int        = base::n_int;
  static constexpr auto       n_frac       = base::n_frac;
  static constexpr auto       n_sign       = base::n_sign;
  static constexpr auto       is_signed    = base::is_signed;
  static constexpr auto       float_factor = base::float_factor;
  static constexpr value_type int_max      = base::int_max;
  static constexpr value_type int_min      = base::int_min;
  static constexpr float_type flt_max      = base::flt_max;
  static constexpr float_type flt_min      = base::flt_min;
  //----------------------------------------------------------------------------
  using base::as_float;
  using base::as_int;
  using base::base;
  using base::fraction;
  using base::load;
  using base::load_float;
  using base::load_int;
  using base::normalize;
  using base::value;
  using base::operator=;
  //----------------------------------------------------------------------------
  template <uint N_intb, uint N_fracb>
  using compatiblefp = fixed_point_int<n_sign, N_intb, N_fracb, true, traits>;
  //----------------------------------------------------------------------------
  template <
    uint N_int_v,
    uint N_frac_v,
    std::enable_if_t<N_int_v <= n_int && N_frac_v <= n_frac>* = nullptr>
  constexpr fixed_point_int (compatiblefp<N_int_v, N_frac_v> fp)
  {
    *this = fp;
  }
  //----------------------------------------------------------------------------
  ~fixed_point_int() = default;
  //----------------------------------------------------------------------------
  // raw constructor
  explicit fixed_point_int (value_type v) { this->_v = v; }
  //----------------------------------------------------------------------------
  // All the operators are lossless
  //----------------------------------------------------------------------------
  template <uint N_sign_v, uint N_int_v, uint N_frac_v>
  constexpr fixed_point_int& operator= (compatiblefp<N_int_v, N_frac_v> fp)
  {
    static_assert (
      N_int_v <= n_int && N_frac_v <= n_frac, "Lossy assignments not allowed");
    *to_base (this) = *to_base (&fp);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator+ (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    constexpr uint op_n_int  = std::max (n_int, N_int_v) + 1;
    constexpr uint op_n_frac = std::max (n_frac, N_frac_v);

    using lossless_type = compatiblefp<op_n_int, op_n_frac>;
    using lvalue_type   = typename lossless_type::value_type;

    constexpr auto shift = (int) (N_frac_v - n_frac);

    lossless_type ret;
    if constexpr (shift > 0) {
      // rhs has more fractional bits
      ret._v = rhs._v;
      ret._v += ashl<shift> ((lvalue_type) this->_v);
    }
    else {
      // instance has more or equal fractional bits
      ret._v = this->_v;
      ret._v += ashl<-shift> ((lvalue_type) rhs._v);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator- (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    constexpr uint op_n_int  = std::max (n_int, N_int_v) + n_sign;
    constexpr uint op_n_frac = std::max (n_frac, N_frac_v);

    using lossless_type = compatiblefp<op_n_int, op_n_frac>;
    using lvalue_type   = typename lossless_type::value_type;

    constexpr auto shift = (int) (N_frac_v - n_frac);

    lossless_type ret;
    if constexpr (shift > 0) {
      // rhs has more fractional bits
      ret._v = ashl<shift> ((lvalue_type) this->_v);
      ret._v -= rhs._v;
    }
    else {
      // instance has more or equal fractional bits
      ret._v = this->_v;
      ret._v -= ashl<-shift> ((lvalue_type) rhs._v);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator* (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    constexpr uint op_n_int  = n_int + N_int_v;
    constexpr uint op_n_frac = n_frac + N_frac_v;

    using lossless_type = compatiblefp<op_n_int, op_n_frac>;
    using lvalue_type   = typename lossless_type::value_type;

    lossless_type ret;
    ret._v = (lvalue_type) this->_v * (lvalue_type) rhs._v;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint N_int_v, uint N_frac_v>
  constexpr auto operator/ (compatiblefp<N_int_v, N_frac_v> rhs) const
  {
    constexpr uint op_n_int  = n_int + N_frac_v;
    constexpr uint op_n_frac = n_frac + N_int_v;

    using lossless_type = compatiblefp<op_n_int, op_n_frac>;

    auto ret = cast<lossless_type>();
    ret._v   = ashl<N_frac_v> (ret._v);
    ret._v /= rhs._v;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T>
  constexpr T cast (T) const
  {
    static_assert (is_fixed_point_int_v<T>);
    static_assert (T::is_signed == is_signed);
    static_assert (std::is_same_v<typename T::promotions, promotions>);

    using ptr       = T*;
    using base_type = std::remove_pointer_t<decltype (to_base (ptr {nullptr}))>;
    T ret;
    ret._v = base::template cast<base_type>()._v;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T>
  constexpr T cast() const
  {
    return cast (T {});
  }
  //----------------------------------------------------------------------------
  static constexpr fixed_point_int from (value_type raw_val)
  {
    return fixed_point_int {raw_val};
  }
  //----------------------------------------------------------------------------
  static constexpr fixed_point_int from_float (float_type flt)
  {
    return fixed_point_int {base::from_float (flt)._v};
  }
  //----------------------------------------------------------------------------
  static constexpr fixed_point_int from_int (value_type intv)
  {
    return fixed_point_int {base::from_int (intv)._v};
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <uint S, uint I, uint F>
  static constexpr auto to_base (fixed_point_int<S, I, F, true, traits>* v)
  {
    return static_cast<detail::fp_int_base<S, I, F>*> (v);
  }
  //----------------------------------------------------------------------------
  template <uint, uint, uint, bool, class>
  friend class fixed_point_int;
};
//------------------------------------------------------------------------------
template <
  class FP,
  uint Src_sign,
  uint Src_int,
  uint Src_frac,
  bool Src_is_lossless,
  class Src_promotions>
constexpr auto fixed_point_int_cast (
  fixed_point_int<Src_sign, Src_int, Src_frac, Src_is_lossless> src)
{
  static_assert (is_fixed_point_int_v<FP>);
  static_assert (std::is_same_v<typename FP::promotions, Src_promotions>);

  using T_src               = typename decltype (src)::value_type;
  using T_dst               = typename FP::value_type;
  constexpr auto bigger_dst = sizeof (T_dst) >= sizeof (T_src);
  using T_big               = std::conditional_t<bigger_dst, T_dst, T_src>;
  constexpr auto shift      = (int) (FP::n_frac - src.frac);

  auto raw = src.template value<T_big>();
  raw      = ash<shift> (raw);
  FP ret;
  ret._v = (T_dst) raw;
  if constexpr (FP::is_lossless && !src.is_lossless) {
    ret.normalize();
  }
  return ret;
}
//------------------------------------------------------------------------------
template <uint S, uint I, uint F, bool L = true, class T = std_fp_types_trait>
using fp_int = fixed_point_int<S, I, F, L, T>;

template <class T, class U>
constexpr auto fp_int_cast (U src)
{
  return fixed_point_int_cast<T> (src);
}
//------------------------------------------------------------------------------
}; // namespace artv
