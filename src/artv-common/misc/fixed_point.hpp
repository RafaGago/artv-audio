#pragma once

#include <algorithm>
#include <assert.h>
#include <gcem.hpp>
#include <type_traits>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec.hpp"

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
    // fp_int_promotion_step<s128, u128, long double>
    >;
  static constexpr uint vector_size = 0;
};
//------------------------------------------------------------------------------
template <uint N>
struct vec_fp_types_trait {
  using promotions = mp_list<
    fp_int_promotion_step<vec<s32, N>, vec<u32, N>, vec<float, N>>,
    fp_int_promotion_step<vec<s64, N>, vec<u64, N>, vec<double, N>>>;
  static constexpr uint vector_size = N;
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
} // namespace detail
//------------------------------------------------------------------------------
template <
  uint N_sign,
  uint N_int,
  uint N_frac,
  bool Is_lossless = true,
  class Traits     = std_fp_types_trait>
class fixpt;

// check if a type is a fixed point int
template <class T>
struct is_fixpt : public std::false_type {};

template <uint N_sign, uint N_int, uint N_frac, bool Is_lossless, class Traits>
struct is_fixpt<fixpt<N_sign, N_int, N_frac, Is_lossless, Traits>>
  : public std::true_type {};

template <class T>
static constexpr bool is_fixpt_v = is_fixpt<T>::value;

template <class T1, class T2>
struct fixpt_are_compatible : public std::false_type {};

// check if two fixed point types are compatible (same sign and promotions)
template <
  uint N_sign1,
  uint N_sign2,
  uint N_int1,
  uint N_int2,
  uint N_frac1,
  uint N_frac2,
  bool Is_lossless1,
  bool Is_lossless2,
  class Traits>
struct fixpt_are_compatible<
  fixpt<N_sign1, N_int1, N_frac1, Is_lossless1, Traits>,
  fixpt<N_sign2, N_int2, N_frac2, Is_lossless2, Traits>>
  : public std::true_type {};

template <class T1, class T2>
static constexpr bool fixpt_are_compatible_v
  = fixpt_are_compatible<T1, T2>::value;

//------------------------------------------------------------------------------
// A fixed point arithmetic class that can use native GCC/clang vector types.
//
// When "Is_lossless" is true the returned types on arithmetic operations expand
// its type to represent the range incresases. In other words, type promotions
// happen.
//
// When "Is_lossless" is false there are no type promotions, but resoultion
// loss.
//
// Notice that for non-lossless types it can make sense to reduce the range as
// much as you know, e.g. using a Q3.7 (16-bit variable), as the (possible)
// range will expand automatically with most arithmetic operations.
//
// "Traits" is a description of the underlying types that can be used, see e.g.
//  "std_fp_types_trait" and "vec_fp_types_trait"
//------------------------------------------------------------------------------
template <uint N_sign, uint N_int, uint N_frac, bool Is_lossless, class Traits>
class fixpt {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_int       = N_int;
  static constexpr uint n_frac      = N_frac;
  static constexpr uint n_sign      = N_sign;
  static constexpr uint n_bits      = n_sign + n_int + n_frac;
  static constexpr uint is_signed   = (n_sign == 1);
  static constexpr bool is_lossless = Is_lossless;
  using traits                      = Traits;
  using promotions                  = typename Traits::promotions;
  //----------------------------------------------------------------------------
private:
  using type_idx
    = mp11::mp_find_if_q<promotions, detail::first_type_n_bits_gt_n<n_bits>>;

  static_assert (n_sign <= 1, "Invalid value. 1 is signed, 0 unsigned");
  static_assert (
    type_idx::value < mp11::mp_size<promotions>::value,
    "Number of bits required not representable by any available type");

  template <class T>
  using enable_if_int_same_signedness = std::enable_if_t<
    std::is_integral_v<T> && std::is_signed_v<T> == (bool) n_sign>;

  template <uint N_intb, uint N_fracb>
  using compatiblefp = fixpt<n_sign, N_intb, N_fracb, is_lossless, traits>;

  using mytype = compatiblefp<n_int, n_frac>;

  template <class T>
  static constexpr bool rhs_is_compatible
    = fixpt_are_compatible_v<T, mytype> && (n_sign == T::n_sign);

  template <class T>
  static constexpr bool rhs_is_assign_compat = rhs_is_compatible<T>
    && (is_lossless || (n_frac >= T::n_frac && n_int >= T::n_int));

  template <class T>
  using enable_if_rhs_is_compat = std::enable_if_t<rhs_is_compatible<T>>;

  template <class T>
  using enable_if_rhs_is_assign_compat
    = std::enable_if_t<rhs_is_assign_compat<T>>;

  template <class T>
  using enable_if_lhs_is_lossless
    = std::enable_if_t<rhs_is_compatible<T> && is_lossless>;

  template <class T>
  using enable_if_lhs_is_lossy
    = std::enable_if_t<rhs_is_compatible<T> && !is_lossless>;

  //----------------------------------------------------------------------------
public:
  using type          = mp11::mp_at<promotions, type_idx>;
  using unsigned_type = typename type::tunsigned;
  using signed_type   = typename type::tsigned;
  using float_type    = typename type::tfloat;
  using value_type = std::conditional_t<is_signed, signed_type, unsigned_type>;

  static constexpr uint vector_size    = traits::vector_size;
  static constexpr uint builtin_sizeof = (vector_size == 0)
    ? sizeof (value_type)
    : sizeof (value_type) / vector_size;

  using builtin_signed   = int_for_size<builtin_sizeof>;
  using builtin_unsigned = uint_for_size<builtin_sizeof>;
  using builtin
    = std::conditional_t<is_signed, builtin_signed, builtin_unsigned>;

  using near_lossless = fixpt<n_sign, n_int, n_frac, true, traits>;
  using near_lossy    = fixpt<n_sign, n_int, n_frac, false, traits>;
  using near_unsigned = fixpt<0, n_int + n_sign, n_frac, is_lossless, traits>;
  using near_signed   = fixpt<
    1,
    !!n_int ? (n_int - n_sign) : 0,
    !n_int ? (n_frac - n_sign) : n_frac,
    is_lossless,
    traits>;
  //----------------------------------------------------------------------------
  constexpr fixpt() { _v = decltype (_v) {}; }

  template <class T, enable_if_rhs_is_assign_compat<T>* = nullptr>
  constexpr fixpt (T other)
  {
    *this = other;
  }
  // TODO, FIX:
  //
  // error: type 'artv::fixpt<1, 0, 15, true> &' cannot be used prior to
  // '::' because it has no members
  //    = fixpt_are_compatible_v<T, mytype> && (n_sign == T::n_sign);
  //
  // fixpt(fixpt&&) = default
  // fixpt& operator= (fixpt&&) = default;
  ~fixpt() = default;
  //----------------------------------------------------------------------------
  // raw value load
  static constexpr fixpt from (value_type raw_val)
  {
    fixpt r;
    r._v = raw_val;
    return r;
  }
  //----------------------------------------------------------------------------
  // lossy conversion from float to fixed point
  static constexpr fixpt from_float (float_type flt)
  {
    // Assertions probably not 100% correct/verified, as the non-constant
    // float/double resolution makes this a non-trivial test. For the
    // meantime this might be enough to catch most mistakes. TODO
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

    fixpt                r;
    float_type           v = flt * float_factor;
    constexpr float_type zero {};
    v += float_type {(v > zero) ? (zero + 0.5f) : (zero - 0.5f)}; // rounding
    r._v = (value_type) v;
    return r;
  }
  //----------------------------------------------------------------------------
  // conversion from integer to fixed point (just the integer part)
  static constexpr fixpt from_int (value_type intv)
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
    return fixpt::from (ashl<n_frac> (intv));
  }
  //----------------------------------------------------------------------------
  constexpr void load (value_type v) { *this = from (v); }
  //----------------------------------------------------------------------------
  constexpr void load_float (float_type v) { *this = from_float (v); }
  //----------------------------------------------------------------------------
  constexpr void load_int (value_type v) { *this = from_int (v); }
  //----------------------------------------------------------------------------
  // cast to any fixp, same caveats that static casting regular integers
  // apply e.g. signedess and range issues.
  template <
    class T,
    std::enable_if_t<fixpt_are_compatible_v<T, mytype>>* = nullptr>
  constexpr T cast (T) const
  {
    // notice that the casts are ignoring the sign...
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
  template <
    class T,
    std::enable_if_t<fixpt_are_compatible_v<T, mytype>>* = nullptr>
  explicit constexpr operator T() const
  {
    return cast<T>();
  }
  //----------------------------------------------------------------------------
  constexpr near_lossless to_lossless() const { return cast<near_lossless>(); };
  //----------------------------------------------------------------------------
  constexpr near_lossy to_lossy() const { return cast<near_lossy>(); };
  //----------------------------------------------------------------------------
  constexpr near_signed to_signed() const { return cast<near_signed>(); };
  //----------------------------------------------------------------------------
  constexpr near_unsigned to_unsigned() const { return cast<near_unsigned>(); };
  //----------------------------------------------------------------------------
  // Add or remove resolution both on the integer and fractional part
  template <int Int, int Frac = 0>
  constexpr auto resize() const
  {
    static_assert (Int >= 0 ? true : (-Int <= n_int));
    static_assert (Frac >= 0 ? true : (-Frac <= n_frac));
    return cast (compatiblefp<n_int + Int, n_frac - Frac> {});
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
  // useful E.g. after random loads to check that the value is in range.
  // Ideally should happen after each load, store, constructor, assignment,
  // etc from external sources, unfortunately this is a numeric/performace
  // sensitive class.
  constexpr void normalize() { _v = alsb_mask (_v, n_frac + n_int); }
  //----------------------------------------------------------------------------
  // raw value
  constexpr value_type value() const { return _v; }
  //----------------------------------------------------------------------------
  template <class T, enable_if_rhs_is_assign_compat<T>* = nullptr>
  constexpr fixpt& operator= (T rhs)
  {
    _v = rhs.cast (*this)._v;
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto& operator+= (T rhs)
  {
    _v = bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a + b; });
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto operator+ (T rhs) const
  {
    auto ret {*this};
    ret += rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossless<T>* = nullptr>
  constexpr auto operator+ (T rhs) const
  {
    constexpr uint int_b  = std::max (n_int, T::n_int) + 1;
    constexpr uint frac_b = std::max (n_frac, T::n_frac);
    return compatiblefp<int_b, frac_b>::from (
      bin_op<int_b, frac_b> (rhs, [] (auto a, auto b) { return a + b; }));
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto& operator-= (T rhs)
  {
    _v = bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a - b; });
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto operator- (T rhs) const
  {
    auto ret {*this};
    ret -= rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossless<T>* = nullptr>
  constexpr auto operator- (T rhs) const
  {
    constexpr uint int_b  = std::max (n_int, T::n_int) + n_sign;
    constexpr uint frac_b = std::max (n_frac, T::n_frac);
    return compatiblefp<int_b, frac_b>::from (
      bin_op<int_b, frac_b> (rhs, [] (auto a, auto b) { return a - b; }));
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto& operator*= (T rhs)
  {
    // only make space to handle (and revert) the shift caused by
    // multiplying.
    constexpr uint int_b  = n_int + T::n_frac;
    constexpr uint frac_b = n_frac;

    using dst_type = compatiblefp<int_b, frac_b>;
    using bigger   = typename dst_type::value_type;
    this->_v
      = (value_type) ashr<T::n_frac> ((bigger) this->_v * (bigger) rhs._v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto operator* (T rhs) const
  {
    auto ret {*this};
    ret *= rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossless<T>* = nullptr>
  constexpr auto operator* (T rhs) const
  {
    constexpr uint int_b  = n_int + T::n_int;
    constexpr uint frac_b = n_frac + T::n_frac;
    using lossless_type   = compatiblefp<int_b, frac_b>;
    using VT              = typename lossless_type::value_type;
    return lossless_type::from ((VT) _v * (VT) rhs._v);
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto& operator/= (T rhs)
  {
    // only make space to handle the pre-shift required for dividing while
    // keeping the same scaling.
    constexpr uint int_b  = n_int + T::n_frac;
    constexpr uint frac_b = n_frac;
    using dst_type        = compatiblefp<int_b, frac_b>;

    auto tmp = cast<dst_type>();
    this->_v = (value_type) (ashl<T::n_frac> (tmp._v) / rhs._v);
    return *this;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossy<T>* = nullptr>
  constexpr auto operator/ (T rhs) const
  {
    auto ret {*this};
    ret /= rhs;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_lhs_is_lossless<T>* = nullptr>
  constexpr auto operator/ (T rhs) const
  {
    constexpr uint int_b  = n_int + T::n_frac;
    constexpr uint frac_b = n_frac + T::n_int;

    using lossless_type = compatiblefp<int_b, frac_b>;

    auto ret = cast<lossless_type>();
    ret._v   = ashl<T::n_frac> (ret._v);
    ret._v /= rhs._v;
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr auto operator-() const
  {
    auto ret {*this};
    ret._v = -_v;
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr auto operator+() const
  {
    auto ret {*this};
    ret._v = +_v;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_rhs_is_compat<T>* = nullptr>
  constexpr auto operator== (T rhs) const
  {
    return bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a == b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_rhs_is_compat<T>* = nullptr>
  constexpr auto operator!= (T rhs) const
  {
    return bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a != b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_rhs_is_compat<T>* = nullptr>
  constexpr auto operator<= (T rhs) const
  {
    return bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a <= b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_rhs_is_compat<T>* = nullptr>
  constexpr auto operator>= (T rhs) const
  {
    return bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a >= b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_rhs_is_compat<T>* = nullptr>
  constexpr auto operator<(T rhs) const
  {
    return bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a < b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_rhs_is_compat<T>* = nullptr>
  constexpr auto operator> (T rhs) const
  {
    return bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a > b; });
  }
  //----------------------------------------------------------------------------
  // Notice: bit shifts are always arithmetic and lossy!
  constexpr fixpt& operator<<= (uint n) { _v = ashl (_v, n); }
  constexpr fixpt& operator>>= (uint n) { _v = ashr (_v, n); }
  constexpr fixpt  operator<< (uint n) { return from (ashl (_v, n)); }
  constexpr fixpt  operator>> (uint n) { return from (ashr (_v, n)); }
  //----------------------------------------------------------------------------
  static constexpr fixpt      max() { return from (raw_max); };
  static constexpr fixpt      min() { return from (raw_min); };
  static constexpr value_type max_raw() { return raw_max; }
  static constexpr value_type min_raw() { return raw_min; }
  static constexpr value_type max_int() { return int_max; }
  static constexpr value_type min_int() { return int_min; }
  static constexpr float_type max_float() { return raw_max; }
  static constexpr float_type min_float() { return raw_min; }
  static constexpr float_type factor() { return float_factor; }
  // TODO: Binary Logic Ops.
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static constexpr unsigned_type unsigned_zeroes = unsigned_type {};
  static constexpr unsigned_type unsigned_ones   = unsigned_type {} - 1U;
  static constexpr unsigned_type unsigned_one    = unsigned_type {} + 1U;

  static constexpr auto float_factor = (float_type) (unsigned_one << n_frac);

  static constexpr value_type raw_max
    = (value_type) (unsigned_ones >> ((builtin_sizeof * 8) - n_bits)) >> n_sign;

  static constexpr value_type raw_min
    = is_signed ? -raw_max - 1 : value_type {};

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
  template <uint Int, uint Frac, class F, class T>
  constexpr auto bin_op (T other, F&& fn) const
  {
    using lossless_type = compatiblefp<Int, Frac>;
    return fn (cast (lossless_type {})._v, other.cast (lossless_type {})._v);
  }
  //----------------------------------------------------------------------------
  template <uint, uint, uint, bool, class>
  friend class fixpt;

  value_type _v;
};
//------------------------------------------------------------------------------
// fixed point type with no promotions
template <uint S, uint I, uint F, class Traits = std_fp_types_trait>
using fixpt_np = fixpt<
  S,
  // round up to the capacity of the builtin type
  (fixpt<S, I, F, false, Traits>::builtin_sizeof * 8) - F - S,
  F,
  false,
  Traits>;

// fixed point type with automatic promotions/range increases
template <uint S, uint I, uint F, class Traits = std_fp_types_trait>
using fixpt_p = fixpt<S, I, F, true, Traits>;
//------------------------------------------------------------------------------
template <bool Sign, class T, class = void>
struct is_fixpt_with_signedness : public std::false_type {};

template <bool Sign, class T>
struct is_fixpt_with_signedness<Sign, T, std::enable_if_t<is_fixpt_v<T>>>
  : public std::integral_constant<bool, T::is_signed == Sign> {};

template <bool Signedness, class T>
static constexpr bool is_fixpt_with_signedness_v
  = is_fixpt_with_signedness<Signedness, T>::value;
//------------------------------------------------------------------------------
// ternary operators can be overloaded. This is mostly done for native vector
// types.
template <class C, class U, std::enable_if_t<is_fixpt_v<U>>* = nullptr>
constexpr auto fixpt_select (C cond, U a, U b)
{
  return U::from (cond ? a.value() : b.value());
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
constexpr auto fixpt_abs (T v)
{
  if constexpr (v.is_signed) {
    constexpr uint shift = (v.builtin_sizeof * 8) - 1;

    auto raw  = v.value();
    auto mask = raw >> shift;
    return T::from ((mask + raw) ^ mask);
  }
  else {
    return v;
  }
}
// TODO: fixpt_ceil, fixpt_floor, fixpt_round, fixpt_exp2...
//------------------------------------------------------------------------------

}; // namespace artv
