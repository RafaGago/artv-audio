#pragma once

#include <algorithm>
#include <cassert>
#include <gcem.hpp>
#include <type_traits>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/comptime_string.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/num.hpp"
#include "artv-common/misc/ratio.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec.hpp"
#include "artv-common/misc/vec_bits.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <class S, class U, class F>
struct fp_int_conversion_step {

  using tsigned   = S;
  using tunsigned = U;
  using tfloat    = F;

  using scalar_signed   = vec_value_type_t<S>;
  using scalar_unsigned = vec_value_type_t<U>;
  using scalar_float    = vec_value_type_t<F>;

  static constexpr size_t n_bits      = sizeof (scalar_signed) * 8;
  static constexpr size_t vector_size = vec_traits_t<tsigned>::size;

  static_assert (sizeof (S) == sizeof (U));
  static_assert (vector_size == vec_traits_t<tfloat>::size);
};
//------------------------------------------------------------------------------
struct std_fp_types_trait {
  using conversions = mp_list<
    fp_int_conversion_step<s8, u8, float>,
    fp_int_conversion_step<s16, u16, float>,
    fp_int_conversion_step<s32, u32, double>,
    fp_int_conversion_step<s64, u64, long double>
    // fp_int_conversion_step<s128, u128, long double>
    >;
};
//------------------------------------------------------------------------------
template <uint N>
struct vec_fp_types_trait {
  using conversions = mp_list<
    fp_int_conversion_step<vec<s32, N>, vec<u32, N>, vec<float, N>>,
    fp_int_conversion_step<vec<s64, N>, vec<u64, N>, vec<double, N>>>;
};
//------------------------------------------------------------------------------
template <class T_signed, class T_unsigned, class T_float, uint Vec_Size = 0>
struct single_fp_type_trait {
  using conversions
    = mp_list<fp_int_conversion_step<T_signed, T_unsigned, T_float>>;
};

//------------------------------------------------------------------------------
namespace detail {
//------------------------------------------------------------------------------
template <uint N>
struct first_type_n_bits_gt_n {
  template <class U>
  using fn = mp11::mp_bool<U::n_bits >= N>;
};

struct fixpt_arith_tag {};

} // namespace detail

//------------------------------------------------------------------------------
// fixed point flags
//------------------------------------------------------------------------------
// dynamic: Arithmetic operations may move the factor (fixed point location)
// and precision. When bit requirement surpasses the capacity of the
// limit/biggest datatype specified in the conversions compilation fails.
//
// They also allow adding signedness as the result of arithmetic operations.
// Non-dynamic types can only operate with types of the same signedness.
static constexpr uint fixpt_dynamic = 1 << 0;
// mixed: Like "dynamic" but tries to truncate all fractional bits before
// failing compilation.
//
// They also allow adding signedness as the result of arithmetic operations, but
// only when bits to allocate are free.
static constexpr uint fixpt_mixed = fixpt_dynamic | (1 << 1);
// rounding: Conversions of fixed point types do round to the nearest instead of
// truncating
static constexpr uint fixpt_rounding = 1 << 2;
// relaxed_frac_assign: Allows assignments where the left hand side would drop
// fractional bits.
static constexpr uint fixpt_relaxed_frac_assign = 1 << 3;
// relaxed_int_assign: Allows assignments where the left hand side would drop
// integer bits.
static constexpr uint fixpt_relaxed_int_assign = 1 << 4;
// allow implicit conversions/arithmetic with floating point types. The floating
// point type will be implicitly converted to the type of the fixed point
// operand in an unsafe fashion
static constexpr uint fixpt_implicit_float = 1 << 5;
// allow implicit conversions/arithmetic with integral point types. The integral
// type will be implicitly converted to the type of the fixed point operand in
// an unsafe fashion.
static constexpr uint fixpt_implicit_int = 1 << 6;
//------------------------------------------------------------------------------
// fixed point flags combinations
//------------------------------------------------------------------------------
// allow implicit conversions/arithmetic with aritmetic types.
static constexpr uint fixpt_implicit
  = fixpt_implicit_float | fixpt_implicit_int;
// relaxed_int_assign: Allows assignments between any fixpt types.
static constexpr uint fixpt_relaxed_assign
  = fixpt_relaxed_frac_assign | fixpt_relaxed_int_assign;
// flexible but dangerous
static constexpr uint fixpt_unsafe = fixpt_implicit | fixpt_relaxed_assign;
//------------------------------------------------------------------------------
template <
  uint N_sign,
  uint N_int,
  uint N_frac,
  uint Flags   = fixpt_dynamic,
  class Traits = std_fp_types_trait>
class fixpt;

// check if a type is a fixed point int
template <class T>
struct is_fixpt : public std::false_type {};

template <uint N_sign, uint N_int, uint N_frac, uint Flags, class Traits>
struct is_fixpt<fixpt<N_sign, N_int, N_frac, Flags, Traits>>
  : public std::true_type {};

template <class T>
static constexpr bool is_fixpt_v = is_fixpt<T>::value;

template <class T1, class T2>
struct fixpt_are_compatible : public std::false_type {};

// check if two fixed point types are compatible (same sign and conversions)
template <
  uint N_sign1,
  uint N_sign2,
  uint N_int1,
  uint N_int2,
  uint N_frac1,
  uint N_frac2,
  uint Flags1,
  uint Flags2,
  class Traits>
struct fixpt_are_compatible<
  fixpt<N_sign1, N_int1, N_frac1, Flags1, Traits>,
  fixpt<N_sign2, N_int2, N_frac2, Flags2, Traits>> : public std::true_type {};

template <class T1, class T2>
static constexpr bool fixpt_are_compatible_v
  = fixpt_are_compatible<T1, T2>::value;

//------------------------------------------------------------------------------
// A fixed point arithmetic class that can use native GCC/clang vector types.
// crafted for audio DSP duties.
//
// See the flags documentation for the capabilities. Also the tests for some
// usage examples.
//
// "Traits" is a compile time-list sorted from smallest to biggest bit
// resolution, enumerating the types that are allowed for doing conversions when
// expanding the type bits or when doing intermediate math and specifying its
// relation to a floating point type. See e.g. "std_fp_types_trait" and
// "vec_fp_types_trait".
//------------------------------------------------------------------------------
template <uint N_sign, uint N_int, uint N_frac, uint Flags, class Traits>
class fixpt {
public:
  //----------------------------------------------------------------------------
  // bit/size related constants
  static constexpr uint n_int     = N_int;
  static constexpr uint n_frac    = N_frac;
  static constexpr uint n_sign    = N_sign;
  static constexpr uint n_bits    = n_sign + n_int + n_frac;
  static constexpr uint is_signed = (n_sign == 1);

  // behavior related constants
  static constexpr uint flags         = Flags;
  static constexpr bool is_saturating = (Flags & fixpt_mixed) != 0;
  static constexpr bool is_dynamic
    = (Flags & fixpt_dynamic) != 0 || is_saturating;
  static constexpr bool rounds_nearest = (Flags & fixpt_rounding) != 0;
  static constexpr bool relaxed_frac_assign
    = (Flags & fixpt_relaxed_frac_assign) != 0;
  static constexpr bool relaxed_int_assign
    = (Flags & fixpt_relaxed_int_assign) != 0;
  static constexpr bool implicit_float = (Flags & fixpt_implicit_float) != 0;
  static constexpr bool implicit_int   = (Flags & fixpt_implicit_int) != 0;

  using traits      = Traits;
  using conversions = typename Traits::conversions;
  //----------------------------------------------------------------------------
private:
  template <class T, uint Fail = -1u>
  static constexpr uint fixpt_sign() noexcept
  {
    if constexpr (is_fixpt_v<T>) {
      return T::n_sign;
    }
    else {
      return Fail;
    }
  }

  template <class T, uint Fail = -1u>
  static constexpr uint fixpt_n_int() noexcept
  {
    if constexpr (is_fixpt_v<T>) {
      return T::n_int;
    }
    else {
      return Fail;
    }
  }

  template <class T, uint Fail = -1u>
  static constexpr uint fixpt_n_frac() noexcept
  {
    if constexpr (is_fixpt_v<T>) {
      return T::n_frac;
    }
    else {
      return Fail;
    }
  }

  template <class T, uint Fail = -1u>
  static constexpr uint fixpt_flags() noexcept
  {
    if constexpr (is_fixpt_v<T>) {
      return T::flags;
    }
    else {
      return Fail;
    }
  }

  static constexpr uint type_idx = mp11::
    mp_find_if_q<conversions, detail::first_type_n_bits_gt_n<n_bits>>::value;
  static constexpr uint n_conversions = mp11::mp_size<conversions>::value;
  static constexpr uint max_conversion_bits
    = mp11::mp_back<conversions>::n_bits;

  static_assert (n_sign <= 1, "Invalid value. 1 is signed, 0 unsigned");
  static_assert (
    type_idx < n_conversions,
    "Number of bits required not representable by any available type");

  using mytype = fixpt<n_sign, n_int, n_frac, flags, traits>;

  // by default unsigned are compatible with signed, but not the other way
  // around
  template <class T>
  static constexpr bool rhs_is_compatible = fixpt_are_compatible_v<T, mytype>;

  template <class T>
  static constexpr bool rhs_is_assign_compat
    = rhs_is_compatible<T> && (n_sign >= fixpt_sign<T>())
    && (relaxed_frac_assign || (n_frac >= fixpt_n_frac<T>()))
    && (relaxed_int_assign || (n_int >= fixpt_n_int<T>()));

  template <class T>
  using enable_if_comparable = std::enable_if_t<rhs_is_compatible<T>>;

  template <class T>
  using enable_if_operand_w_dynamic
    = std::enable_if_t<rhs_is_compatible<T> && is_dynamic>;

  template <class T>
  using enable_if_operand_w_static = std::enable_if_t<
    rhs_is_compatible<T> && (n_sign >= fixpt_sign<T>()) && !is_dynamic>;
  //----------------------------------------------------------------------------
  template <uint N_intb, uint N_fracb, uint N_signv, uint Flagsv>
  static constexpr auto convert_impl()
  {
    constexpr bool saturating = (Flagsv & fixpt_mixed) != 0;
    constexpr bool dynamic    = (Flagsv & fixpt_dynamic) != 0 || is_saturating;

    if constexpr (dynamic && saturating) {
      return fixpt<
        N_signv,
        N_intb, // the most significant integer bits are always respected
        std::min (N_fracb, max_conversion_bits - N_signv - N_intb),
        Flagsv,
        traits> {};
    }
    else if constexpr (is_dynamic) {
      return fixpt<N_signv, N_intb, N_fracb, Flagsv, traits> {};
    }
    else {
      // static
      using cdynamic
        = fixpt<N_signv, N_intb, N_fracb, Flagsv | fixpt_dynamic, traits>;
      return fixpt<
        N_signv,
        cdynamic::scalar_bits - N_signv - N_fracb,
        N_fracb,
        Flagsv,
        traits> {};
    }
  }
  //----------------------------------------------------------------------------
  template <uint C_Flags>
  using convert_flags = fixpt<n_sign, n_int, n_frac, C_Flags, traits>;

  template <
    uint N_intb,
    uint N_fracb,
    uint N_signv = n_sign,
    uint Flagsv  = flags>
  using convert_raw
    = decltype (convert_impl<N_intb, N_fracb, N_signv, Flagsv>());
  //----------------------------------------------------------------------------
public:
  using type = mp11::mp_at_c<conversions, type_idx>;
  // types next might or not be vectors
  using uint_type  = typename type::tunsigned;
  using sint_type  = typename type::tsigned;
  using float_type = typename type::tfloat;
  // types next are obviously scalars
  using scalar_sint  = typename type::scalar_signed;
  using scalar_uint  = typename type::scalar_unsigned;
  using scalar_float = typename type::scalar_float;

  using value_type  = std::conditional_t<is_signed, sint_type, uint_type>;
  using scalar_type = std::conditional_t<is_signed, scalar_sint, scalar_uint>;
  // vector_size == 0 -> scalar
  static constexpr uint vector_size = type::vector_size;
  static constexpr uint scalar_bits = type::n_bits;

  static_assert (
    is_dynamic || n_bits == scalar_bits,
    "Static/lossy types have to use all bits of the underlying integral type");

  // convert to a type with a different amount of bits but the same signedness
  // and behavior
  template <uint N_intb, uint N_fracb, uint N_signv = n_sign>
  using convert = convert_raw<N_intb, N_fracb, N_signv>;

  // sign related instatiations
  using unsigned_twin = fixpt<0, n_int + n_sign, n_frac, flags, traits>;
  using signed_twin   = fixpt<
    1,
    !!n_int ? (n_int - n_sign) : 0,
    !n_int ? (n_frac - n_sign) : n_frac,
    flags,
    traits>;
  using with_sign = fixpt<1, n_int, n_frac, flags, traits>;

  // other fixpt instantiations with the same bit flags but with different
  // behavior flags
  using dynamic_twin      = convert_flags<(flags & ~uint {3}) | fixpt_dynamic>;
  using static_twin       = convert_flags<flags & ~uint {3}>;
  using mixed_twin        = convert_flags<(flags & ~uint {3}) | fixpt_mixed>;
  using rounding_twin     = convert_flags<flags | fixpt_rounding>;
  using truncating_twin   = convert_flags<flags & ~fixpt_rounding>;
  using strict_frac_twin  = convert_flags<flags & ~fixpt_relaxed_frac_assign>;
  using relaxed_frac_twin = convert_flags<flags | fixpt_relaxed_frac_assign>;
  using strict_int_twin   = convert_flags<flags & ~fixpt_relaxed_int_assign>;
  using relaxed_int_twin  = convert_flags<flags | fixpt_relaxed_int_assign>;
  using strict_twin       = convert_flags<flags & ~fixpt_relaxed_assign>;
  using relaxed_twin      = convert_flags<flags | fixpt_relaxed_assign>;
  using implicit_float_twin = convert_flags<flags | fixpt_implicit_float>;
  using no_float_twin       = convert_flags<flags & ~fixpt_implicit_float>;
  using implicit_int_twin   = convert_flags<flags | fixpt_implicit_int>;
  using no_int_twin         = convert_flags<flags & ~fixpt_implicit_int>;
  using implicit_twin       = convert_flags<flags | fixpt_implicit>;
  using no_implicit_twin    = convert_flags<flags & ~fixpt_implicit>;
  //----------------------------------------------------------------------------
  constexpr fixpt() noexcept { _v = decltype (_v) {}; }

  template <class T, std::enable_if_t<rhs_is_assign_compat<T>>* = nullptr>
  constexpr fixpt (T other) noexcept
  {
    *this = other;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_vec_or_scalar_t<T>* = nullptr>
  constexpr fixpt (num<T> v) noexcept
  {
    *this = v;
  }
  //----------------------------------------------------------------------------
  // implicit int conversion, notice that it might fail on the operator= impl
  template <class T, enable_if_vec_or_scalar_t<T>* = nullptr>
  constexpr fixpt (T v, detail::fixpt_arith_tag) noexcept
  {
    *this = v;
  }
  //----------------------------------------------------------------------------
  // explicit raw load of the value type
  constexpr explicit fixpt (value_type v) noexcept { _v = v; }
  //----------------------------------------------------------------------------
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
  static constexpr fixpt from (value_type raw_val) noexcept
  {
    return fixpt {raw_val};
  }
  //----------------------------------------------------------------------------
  // lossy conversion from float to fixed point
  template <
    class T,
    std::enable_if_t<
      std::is_same_v<T, float_type> || std::is_floating_point_v<T>>* = nullptr>
  static constexpr fixpt from_float (T flt) noexcept
  {
    // Assertions probably not 100% correct/verified, as the non-constant
    // float/double resolution makes this a non-trivial test. For the
    // meantime this might be enough to catch most mistakes. TODO

    auto mmax = max_float();
    auto mmin = min_float();

    auto lte_max = (flt <= max_float());
    auto gte_min = (flt >= min_float());

    if constexpr (vector_size == 0) {
      assert (lte_max);
      assert (gte_min);
    }
    else {
      for (uint i = 0; i < vector_size; ++i) {
        assert (lte_max[i]);
        assert (gte_min[i]);
      }
    }

    fixpt                r;
    float_type           v = vec_cast<scalar_float> (flt * factor());
    constexpr float_type zero {};
    v += float_type {(v > zero) ? (zero + 0.5f) : (zero - 0.5f)}; // rounding
    r._v = (value_type) v;
    return r;
  }
  //----------------------------------------------------------------------------
  // conversion from integer to fixed point (just the integer part)
  template <
    class T,
    std::enable_if_t<
      std::is_same_v<T, value_type> || std::is_integral_v<T>>* = nullptr>
  static constexpr fixpt from_int (T intv) noexcept
  {
    auto casted = vec_set<vector_size> (vec_cast<scalar_type> (intv));
    if constexpr (vector_size == 0) {
      assert (casted <= max_int());
      assert (casted >= min_int());
      assert (intv == (T) casted);
    }
    else {
      assert (vec_is_all_ones (casted <= max_int()));
      assert (vec_is_all_ones (casted >= min_int()));
      assert (vec_is_all_ones (intv == vec_cast<T> (casted)));
    }
    return fixpt::from (ashl<n_frac> (casted));
  }
  //----------------------------------------------------------------------------
  constexpr void load (value_type v) noexcept { *this = from (v); }
  //----------------------------------------------------------------------------
  constexpr void load_float (float_type v) noexcept { *this = from_float (v); }
  //----------------------------------------------------------------------------
  constexpr void load_int (value_type v) noexcept { *this = from_int (v); }
  //----------------------------------------------------------------------------
  // cast to any fixp, same caveats that static casting regular integers
  // apply e.g. signedness and range issues might not be respected.
  template <
    class T,
    std::enable_if_t<fixpt_are_compatible_v<T, mytype>>* = nullptr>
  constexpr T cast (T) const noexcept
  {
    // notice that the casts are ignoring the sign...
    constexpr auto shift = (int) (T::n_frac - n_frac);
    T              ret;
    using scalar = typename T::scalar_type;
    if constexpr (shift > 0) {
      // dst has more fractional bits
      ret._v = ashl<shift> (vec_cast<scalar> (_v));
    }
    else if constexpr (shift < 0) {
      ret._v = vec_cast<scalar> (
        ret.template ashr_truncate<-shift, scalar_type> (_v));
    }
    else {
      ret._v = vec_cast<scalar> (_v);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  template <
    class T,
    std::enable_if_t<fixpt_are_compatible_v<T, mytype>>* = nullptr>
  constexpr T cast() const noexcept
  {
    return cast (T {});
  }
  //----------------------------------------------------------------------------
  // cast by capturing the specificatio of T, but keeping the current traits
  // (promotions/conversion types), as traits conversion is not possible
  template <class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  constexpr auto spec_cast (T) const noexcept
  {
    return cast (fixpt<T::n_sign, T::n_int, T::n_frac, T::flags, traits> {});
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  constexpr auto spec_cast() const noexcept
  {
    return spec_cast (T {});
  }
  //----------------------------------------------------------------------------
  // casting flags may result in invalid sizes for static types (they might
  // expand the integer part).
  template <uint Flags_v>
  constexpr auto flags_cast() const noexcept
  {
    return cast (convert_raw<n_int, n_frac, n_sign, Flags_v> {});
  }
  //----------------------------------------------------------------------------
  // Shorthand functions to just convert flags (could be considered bloat), Just
  // for being able to write better looking client code.
  //----------------------------------------------------------------------------
  constexpr dynamic_twin to_dynamic() const noexcept
  {
    return flags_cast<dynamic_twin::flags>();
  };
  constexpr static_twin to_static() const noexcept
  {
    return flags_cast<static_twin::flags>();
  };
  constexpr mixed_twin to_mixed() const noexcept
  {
    return flags_cast<mixed_twin::flags>();
  };
  //----------------------------------------------------------------------------
  constexpr truncating_twin to_truncating() const noexcept
  {
    return cast<truncating_twin>();
  };
  constexpr rounding_twin to_rounding() const noexcept
  {
    return cast<rounding_twin>();
  };
  //----------------------------------------------------------------------------
  constexpr relaxed_frac_twin to_relaxed_frac_assign() const noexcept
  {
    return cast<relaxed_frac_twin>();
  }
  constexpr strict_frac_twin to_strict_frac_assign() const noexcept
  {
    return cast<strict_frac_twin>();
  }
  constexpr relaxed_int_twin to_relaxed_int_assign() const noexcept
  {
    return cast<relaxed_int_twin>();
  }
  constexpr strict_int_twin to_strict_int_assign() const noexcept
  {
    return cast<strict_int_twin>();
  }
  constexpr relaxed_twin to_relaxed_assign() const noexcept
  {
    return cast<relaxed_twin>();
  }
  constexpr strict_twin to_strict_assign() const noexcept
  {
    return cast<strict_twin>();
  }
  //----------------------------------------------------------------------------
  constexpr implicit_float_twin to_implicit_float() const noexcept
  {
    return cast<implicit_float_twin>;
  }
  constexpr no_float_twin to_no_float() const noexcept
  {
    return cast<no_float_twin>;
  }
  constexpr implicit_int_twin to_implicit_int() const noexcept
  {
    return cast<implicit_int_twin>;
  }
  constexpr no_int_twin to_no_int() const noexcept { return cast<no_int_twin>; }

  constexpr implicit_twin to_implicit_arith() const noexcept
  {
    return cast<implicit_twin>;
  }
  constexpr no_implicit_twin to_no_implicit_arith() const noexcept
  {
    return cast<no_implicit_twin>;
  }
  //----------------------------------------------------------------------------
  constexpr signed_twin to_signed() const noexcept
  {
    return cast<signed_twin>();
  };
  constexpr unsigned_twin to_unsigned() const noexcept
  {
    return cast<unsigned_twin>();
  };
  constexpr with_sign add_sign() const noexcept { return cast<with_sign>(); };
  //----------------------------------------------------------------------------
  // Add or remove resolution both on the integer and fractional part. Negative
  // values remove resolution.
  template <int Int, int Frac = 0>
  constexpr auto resize() const noexcept
  {
    static_assert (Int >= 0 ? true : (-Int <= n_int));
    static_assert (Frac >= 0 ? true : (-Frac <= n_frac));
    return cast (convert<n_int + Int, n_frac + Frac> {});
  }
  // Set resolution in absolute terms
  template <uint Int, uint Frac>
  constexpr auto set_size() const noexcept
  {
    return cast (convert<Int, Frac> {});
  }
  //----------------------------------------------------------------------------
  // Get value as a float
  constexpr float_type to_floatp() const noexcept
  {
    constexpr scalar_float factor = (scalar_float) 1. / float_factor;
    return factor * vec_cast<scalar_float> (_v);
  }
  //----------------------------------------------------------------------------
  // Get value as a an int (equivalent to floor)
  constexpr value_type to_int() const noexcept { return ashr<n_frac> (_v); }
  //----------------------------------------------------------------------------
  constexpr value_type round()
  {
    return to_rounding().template set_size<n_int, 0>().to_int();
  }
  //----------------------------------------------------------------------------
  // returns the fractional part only. for signed types it might be negative.
  constexpr fixpt<n_sign, 0, n_frac, flags, traits> fractional() const noexcept
  {
    auto v = _v;
    if constexpr (is_signed) {
      if (v >= 0) {
        v &= frac_mask;
      }
      else {
        v |= ~frac_mask;
      }
    }
    else {
      v &= frac_mask;
    }
    return fixpt<n_sign, 0, n_frac, flags, traits>::from (v);
  }
  //----------------------------------------------------------------------------
  // useful E.g. after random loads to check that the value is in range.
  // Ideally should happen after each load, store, constructor, assignment,
  // etc from external sources, unfortunately this is a numeric/performace
  // sensitive class.
  constexpr void normalize() noexcept { _v = alsb_mask (_v, n_bits); }
  //----------------------------------------------------------------------------
  constexpr bool is_normalized() const noexcept
  {
    auto v = (_v >= raw_min) && (_v <= raw_max);
    if constexpr (vector_size == 0) {
      return v;
    }
    else {
      return vec_is_all_ones (v);
    }
  }
  //----------------------------------------------------------------------------
  // complement-2 signed integers have more range on the negative part by 1 ulp.
  // This function clamps the negative part
  constexpr fixpt symmetric_clamp() noexcept
  {
    if constexpr (is_signed) {
      return from (_v ^ ((_v == raw_min) & 1));
    }
    else {
      return *this;
    }
  }
  //----------------------------------------------------------------------------
  // raw underlying integer value(s)
  constexpr value_type value() const noexcept { return _v; }
  //----------------------------------------------------------------------------
  // maximum value representable
  static constexpr fixpt max() noexcept
  {
    return from (vec_set<vector_size> (raw_max));
  };
  // minimum value representable
  static constexpr fixpt min() noexcept
  {
    return from (vec_set<vector_size> (raw_min));
  };
  // value of the step. basically 1 on the underlying integer type.
  static constexpr fixpt epsilon() noexcept
  {
    return from (vec_set<vector_size> ((scalar_type) 1));
  }
  // maximum value representable as a raw integer representation. Equivalent to
  // max().value()
  static constexpr value_type max_raw() noexcept
  {
    return vec_set<vector_size> (raw_max);
  }
  // maximum value representable as a raw integer representation. Equivalent to
  // min().value()
  static constexpr value_type min_raw() noexcept
  {
    return vec_set<vector_size> (raw_min);
  }
  // floor of "max_raw()" on the internal representation (zeroes at the right)
  static constexpr value_type max_int() noexcept
  {
    return vec_set<vector_size> (int_max);
  }
  // floor of "min_raw()" on the internal representation (zeroes at the right)
  static constexpr value_type min_int() noexcept
  {
    return vec_set<vector_size> (int_min);
  }
  // Maximum value representable as float. Beware that depending on the traits
  // configuration there might be resolution issues because of the provided
  // floating point type may not have enough mantissa bits to represent it.
  static constexpr float_type max_float() noexcept
  {
    return vec_set<vector_size> (flt_max);
  }
  // see "max_float"
  static constexpr float_type min_float() noexcept
  {
    return vec_set<vector_size> (flt_min);
  }
  // The location of the fixed point normalized. In other words, the factor that
  // makes the raw integer representation to result in the floating point value.
  static constexpr float_type factor() noexcept
  {
    return vec_set<vector_size> (float_factor);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<rhs_is_assign_compat<T>>* = nullptr>
  constexpr fixpt& operator= (T rhs) noexcept
  {
    _v = rhs.cast (*this)._v;
    assert (
      to_int() == rhs.to_int()
      || ((relaxed_frac_assign || relaxed_int_assign) && rounds_nearest));
    return *this;
  }

  template <
    class T,
    std::enable_if_t<is_fixpt_v<T> && !rhs_is_assign_compat<T>>* = nullptr>
  constexpr fixpt& operator= (T rhs) noexcept
  {
    if constexpr (std::is_same_v<typename T::traits, traits>) {
      static_assert (
        sizeof (T) == 0,
        "Assignment loses range or resolution, manually cast if suitable");
    }
    else {
      static_assert (
        sizeof (T) == 0,
        "Can't assign to a fixed point type with different conversions");
    }
  }

  template <class T, enable_if_any_int_vec_or_scalar_t<T>* = nullptr>
  constexpr fixpt& operator= (num<T> rhs) noexcept
  {
    *this = from_int (rhs.value);
    return *this;
  }

  template <class T, enable_if_floatpt_vec_or_scalar_t<T>* = nullptr>
  constexpr fixpt& operator= (num<T> rhs) noexcept
  {
    *this = from_float (rhs.value);
    return *this;
  }

  template <class T, enable_if_any_int_vec_or_scalar_t<T>* = nullptr>
  constexpr fixpt& operator= (T rhs) noexcept
  {
    if constexpr (implicit_int) {
      *this = from_int (rhs);
      return *this;
    }
    else {
      static_assert (
        sizeof (T) == 0,
        "Implicit assignments/conversions to integral types not allowed");
      return *this;
    }
  }

  template <class T, enable_if_floatpt_vec_or_scalar_t<T>* = nullptr>
  constexpr fixpt& operator= (T rhs) noexcept
  {
    if constexpr (implicit_float) {
      *this = from_float (rhs);
      return *this;
    }
    else {
      static_assert (
        sizeof (T) == 0,
        "Implicit assignments/conversions to floating point types not allowed");
      return *this;
    }
  }
  //----------------------------------------------------------------------------
  template <
    class T,
    std::enable_if_t<fixpt_are_compatible_v<T, mytype>>* = nullptr>
  explicit constexpr operator T() const noexcept
  {
    return cast (T {});
  }
  //----------------------------------------------------------------------------
  template <
    class T,
    std::enable_if_t<
      std::is_same_v<T, float_type> || std::is_floating_point_v<T>>* = nullptr>
  explicit constexpr operator T() const noexcept
  {
    return vec_cast<vec_value_type_t<T>> (to_floatp());
  }
  //----------------------------------------------------------------------------
  template <
    class T,
    std::enable_if_t<
      std::is_same_v<T, value_type> || std::is_integral_v<T>>* = nullptr>
  explicit constexpr operator T() const noexcept
  {
    return vec_cast<vec_value_type_t<T>> (to_int());
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto& operator+= (T rhs) noexcept
  {
    _v = bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a + b; });
    return *this;
  }

  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto operator+ (T rhs) const noexcept
  {
    auto ret {*this};
    ret += rhs;
    return ret;
  }

  template <class T, enable_if_operand_w_dynamic<T>* = nullptr>
  constexpr auto operator+ (T rhs) const noexcept
  {
    assert (rhs.is_normalized()); // making sure it is in range
    constexpr uint int_b  = std::max (n_int, T::n_int) + 1;
    constexpr uint frac_b = std::max (n_frac, T::n_frac);
    constexpr uint sign   = n_sign | T::n_sign;

    (void) get_and_assert_saturation<int_b, frac_b, sign>();

    return convert<int_b, frac_b, sign>::from (
      bin_op<int_b, frac_b, sign> (rhs, [] (auto a, auto b) { return a + b; }));
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto& operator-= (T rhs) noexcept
  {
    _v = bin_op<n_int, n_frac> (rhs, [] (auto a, auto b) { return a - b; });
    return *this;
  }

  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto operator- (T rhs) const noexcept
  {
    auto ret {*this};
    ret -= rhs;
    return ret;
  }

  template <class T, enable_if_operand_w_dynamic<T>* = nullptr>
  constexpr auto operator- (T rhs) const noexcept
  {
    assert (rhs.is_normalized()); // making sure it is in range
    constexpr uint sign   = n_sign | T::n_sign;
    constexpr uint int_b  = std::max (n_int, T::n_int) + sign;
    constexpr uint frac_b = std::max (n_frac, T::n_frac);

    (void) get_and_assert_saturation<int_b, frac_b, sign>();

    return convert<int_b, frac_b, sign>::from (
      bin_op<int_b, frac_b, sign> (rhs, [] (auto a, auto b) { return a - b; }));
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto& operator*= (T rhs) noexcept
  {
    // only make space to handle (and revert) the shift caused by
    // multiplying.
    constexpr uint int_b  = n_int + T::n_frac;
    constexpr uint frac_b = n_frac;

    using converted  = convert<int_b, frac_b>;
    using new_scalar = typename converted::scalar_type;
    auto mul = vec_cast<new_scalar> (_v) * vec_cast<new_scalar> (rhs._v);

    _v = vec_cast<scalar_type> (ashr_truncate<T::n_frac, new_scalar> (mul));
    return *this;
  }

  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto operator* (T rhs) const noexcept
  {
    auto ret {*this};
    ret *= rhs;
    return ret;
  }

  template <class T, enable_if_operand_w_dynamic<T>* = nullptr>
  constexpr auto operator* (T rhs) const noexcept
  {
    assert (rhs.is_normalized()); // making sure it is in range

    constexpr uint int_b  = n_int + T::n_int;
    constexpr uint frac_b = n_frac + T::n_frac;
    constexpr uint sign   = n_sign | T::n_sign;

    constexpr auto n_sat = get_and_assert_saturation<int_b, frac_b, sign>();

    using converted = convert<int_b, frac_b, sign>;
    using scalar    = typename converted::scalar_type;
    if constexpr (n_sat > 0) {
      // For mixed types: dropping fractional bits before considering that not
      // enough datatype resolution is available. Notice that the
      // "get_and_assert_saturation" will "static_assert" when there are no more
      // bits available.
      constexpr uint n_frac_kept = converted::n_frac;
      // try with a balanced amount of fractional bits.
      constexpr uint kept1 = (n_frac_kept + 1) / 2; // div_ceil
      constexpr uint kept2 = n_frac_kept / 2;

      constexpr uint kept_lhs = (n_frac > T::n_frac) ? kept1 : kept2;
      constexpr uint kept_rhs = (n_frac > T::n_frac) ? kept2 : kept1;

      if constexpr ((n_frac >= kept_lhs) && (T::n_frac >= kept_rhs)) {
        // balanced fractional bit dropping
        constexpr uint lhs_drop = n_frac - kept_lhs;
        constexpr uint rhs_drop = T::n_frac - kept_rhs;

        return converted::from (
          vec_cast<scalar> (ashr_truncate_ns<lhs_drop> (_v))
          * vec_cast<scalar> (ashr_truncate_ns<rhs_drop> (rhs._v)));
      }
      else if constexpr ((n_frac >= kept_lhs)) {
        // drop more bits from lhs.
        constexpr uint kept_rhs_debt = kept_rhs - T::n_frac;
        constexpr uint lhs_drop      = n_frac - kept_lhs - kept_rhs_debt;

        return converted::from (
          vec_cast<scalar> (ashr_truncate_ns<lhs_drop> (_v))
          * vec_cast<scalar> (rhs._v));
      }
      else if constexpr ((n_frac >= kept_lhs) && (T::n_frac >= kept_rhs)) {
        // drop more bits from rhs.
        constexpr uint kept_lhs_debt = kept_lhs - n_frac;
        constexpr uint rhs_drop      = T::n_frac - kept_rhs - kept_lhs_debt;

        return converted::from (
          vec_cast<scalar> (_v)
          * vec_cast<scalar> (ashr_truncate_ns<rhs_drop> (rhs._v)));
      }
      else {
        static_assert (sizeof (T) == 0, "Unreachable!");
      }
    }
    else {
      return converted::from (
        vec_cast<scalar> (_v) * vec_cast<scalar> (rhs._v));
    }
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto& operator/= (T rhs) noexcept
  {
    // only make space to handle the pre-shift required for dividing while
    // keeping the same scaling.
    constexpr uint int_b  = n_int + T::n_frac;
    constexpr uint frac_b = n_frac;

    using converted = convert<int_b, frac_b>;

    auto tmp = cast<converted>();
    this->_v = vec_cast<scalar_type> (ashl<T::n_frac> (tmp._v) / rhs._v);
    return *this;
  }

  template <class T, enable_if_operand_w_static<T>* = nullptr>
  constexpr auto operator/ (T rhs) const noexcept
  {
    auto ret {*this};
    ret /= rhs;
    return ret;
  }

  template <class T, enable_if_operand_w_dynamic<T>* = nullptr>
  constexpr auto operator/ (T rhs) const noexcept
  {
    // reminder: The easiest explanation for me of fixed point division is based
    // on considering the fact that a division by a decimal number is a
    // multiplication. E.g
    //
    // 1 / 0.001 = 1000;
    //
    // Considering the properties of the integer value 1 in a fixed point int:
    //
    // - On fixed point it is the smallest (positive) number representable.
    // - On fixed point, when used as a divisor it has to multiply the
    //   numerator/dividend.
    // - On integers, dividing an integer by 1 keeps the result is the same.
    //
    // So it becomes obvious that fixed point division requires adjusting the
    // numerator/dividend in a way such as when it is divided by one it is
    // equivalent to a multiplication. To do so it has to be shifted left from
    // the number of fractional bits of the divisor relative to the final
    // desired fixed point location.

    assert (rhs.is_normalized()); // make sure it is in range

    constexpr uint int_b  = n_int + T::n_frac;
    constexpr uint frac_b = n_frac + T::n_int;
    constexpr uint sign   = n_sign | T::n_sign;

    using converted = convert<int_b, frac_b, sign>;
    using scalar    = typename converted::scalar_type;

    constexpr int  prediv_req_pos = converted::n_frac + T::n_frac;
    constexpr auto curr_pos       = (int) n_frac;
    constexpr auto shift          = prediv_req_pos - curr_pos;

    get_and_assert_saturation<int_b, frac_b, sign>();

    auto v = vec_cast<scalar> (_v);
    v      = ::artv::ash<shift> (v);
    return converted {v / rhs._v};
  }
  //----------------------------------------------------------------------------
  constexpr auto operator-() const noexcept
  {
    auto ret {*this};
    ret._v = -_v;
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr auto operator+() const noexcept
  {
    auto ret {*this};
    ret._v = +_v;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_comparable<T>* = nullptr>
  constexpr auto operator== (T rhs) const noexcept
  {
    return cmp_op (rhs, [] (auto a, auto b) { return a == b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_comparable<T>* = nullptr>
  constexpr auto operator!= (T rhs) const noexcept
  {
    return cmp_op (rhs, [] (auto a, auto b) { return a != b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_comparable<T>* = nullptr>
  constexpr auto operator<= (T rhs) const noexcept
  {
    return cmp_op (rhs, [] (auto a, auto b) { return a <= b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_comparable<T>* = nullptr>
  constexpr auto operator>= (T rhs) const noexcept
  {
    return cmp_op (rhs, [] (auto a, auto b) { return a >= b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_comparable<T>* = nullptr>
  constexpr auto operator<(T rhs) const noexcept
  {
    return cmp_op (rhs, [] (auto a, auto b) { return a < b; });
  }
  //----------------------------------------------------------------------------
  template <class T, enable_if_comparable<T>* = nullptr>
  constexpr auto operator> (T rhs) const noexcept
  {
    return cmp_op (rhs, [] (auto a, auto b) { return a > b; });
  }
  //----------------------------------------------------------------------------
  // Notice: bit shifts are always arithmetic and lossy!
  constexpr fixpt& operator<<= (uint n) noexcept { _v = ashl (_v, n); }
  constexpr fixpt& operator>>= (uint n) noexcept { _v = ashr (_v, n); }
  constexpr fixpt  operator<< (uint n) noexcept { return from (ashl (_v, n)); }
  constexpr fixpt  operator>> (uint n) noexcept { return from (ashr (_v, n)); }
  //----------------------------------------------------------------------------
  // arithmetic shift via floating point movement. positive values of N shift
  // left. Dynamic type aware.
  template <int N>
  constexpr auto ash()
  {
    if constexpr (is_dynamic) {
      if constexpr (N >= 0) {
        constexpr uint new_int = n_int + N;
        if constexpr (n_frac >= N) {
          // no-op,.just a virtual point change
          return convert<new_int, n_frac - N> {_v};
        }
        else {
          // requires adding some zeros at the right
          auto           big       = convert<new_int, 0> {_v};
          constexpr uint shift_rem = N - n_frac;
          return big << shift_rem;
        }
      }
      else { // N is negative, shift right
        constexpr int new_frac = (int) n_frac - N;
        constexpr int new_int  = (n_int >= -N) ? (int) n_int + N : 0;
        // no-op,.just a virtual point change (or adding virtual zeroes at the
        // left)
        return convert<new_int, new_frac> {_v};
      }
    }
    else {
      auto r = *this;
      r._v   = ash<N> (r._v);
      return r;
    }
  }
  //----------------------------------------------------------------------------
  // TODO: Binary Logic Ops. Do they make sense.
  //----------------------------------------------------------------------------
private:
  template <uint G_int, uint G_frac, uint Signv>
  static constexpr uint get_and_assert_saturation()
  {
    constexpr uint required_bits = Signv + G_int + G_frac;
    if constexpr (is_saturating) {
      static_assert (
        (required_bits - G_frac) <= max_conversion_bits,
        "Arithmetic operation needs more integer bits than the biggest convertible type supports");
      return (required_bits > max_conversion_bits)
        ? required_bits - max_conversion_bits
        : 0;
    }
    else if constexpr (is_dynamic && (required_bits > max_conversion_bits)) {
      static_assert (
        (required_bits - G_frac) <= max_conversion_bits,
        "Arithmetic operation needs more integer bits than the biggest convertible type supports");
      static_assert (
        (required_bits - G_frac) > max_conversion_bits,
        "Arithmetic operation needs more fractional bits than the biggest convertible type supports");
    }
    return 0;
  }
  //----------------------------------------------------------------------------
  template <uint N_drop, class T_scalar, class T>
  static constexpr auto ashr_truncate (T v)
  {
    if constexpr (rounds_nearest && (N_drop > 0)) {
      constexpr auto roundval = bit<T_scalar> (N_drop - 1);
      if constexpr (std::is_signed_v<T_scalar>) {
        v += (v >= 0) ? roundval : -roundval;
      }
      else {
        v += roundval;
      }
    }
    return ashr<N_drop> (v);
  }
  //----------------------------------------------------------------------------
  template <uint N_drop, class T>
  static constexpr auto ashr_truncate_ns (T v)
  {
    return ashr_truncate<N_drop, T> (v);
  }
  //----------------------------------------------------------------------------
  static constexpr auto frac_mask
    = n_frac != 0 ? lsb_mask<scalar_uint> (n_frac) : scalar_uint {};

  static constexpr scalar_uint unsigned_zeroes = ((scalar_uint) 0);
  static constexpr scalar_uint unsigned_ones   = ((scalar_uint) 0) - 1u;
  static constexpr scalar_uint unsigned_one    = ((scalar_uint) 0) + 1u;

  static constexpr scalar_type raw_max
    = (scalar_type) ((unsigned_ones >> (scalar_bits - n_bits)) >> n_sign);

  static constexpr scalar_type raw_min
    = is_signed ? -raw_max - 1 : scalar_type {};

  static constexpr scalar_type int_max = (n_frac >= scalar_bits)
    ? scalar_type {}
    : (scalar_type) (((scalar_uint) raw_max) >> n_frac);

  static constexpr scalar_type int_min
    = is_signed ? -int_max - 1 : scalar_type {};

  // these limits might not always be precise depending on the floating point
  // type
  static constexpr auto float_factor = (n_frac < scalar_bits)
    ? (scalar_float) (unsigned_one << n_frac)
    : ((scalar_float) int_max) + scalar_float {1};

  static constexpr scalar_float flt_max = (scalar_float) raw_max / float_factor;
  static constexpr scalar_float flt_min = (scalar_float) raw_min / float_factor;
  //----------------------------------------------------------------------------
  template <uint Int, uint Frac, uint Signv, class F, class T>
  constexpr auto bin_op (T other, F&& fn) const noexcept
  {
    using converted = convert<Int, Frac, Signv>;
    return fn (cast (converted {})._v, other.cast (converted {})._v);
  }
  //----------------------------------------------------------------------------
  template <uint Int, uint Frac, class F, class T>
  constexpr auto bin_op (T other, F&& fn) const noexcept
  {
    return bin_op<Int, Frac, n_sign> (other, std::forward<F> (fn));
  }
  //----------------------------------------------------------------------------
  template <class T, class F>
  constexpr auto cmp_op (T other, F&& fn) const noexcept
  {
    static_assert (
      n_sign == T::n_sign,
      "Signedness mismatch. Manually cast before comparison.");
    constexpr uint intb  = std::max (n_int, T::n_int);
    constexpr uint fracb = std::max (n_frac, T::n_frac);
    using converted      = convert<intb, fracb>;
    return fn (cast (converted {})._v, other.cast (converted {})._v);
  }
  //----------------------------------------------------------------------------
  template <uint, uint, uint, uint, class>
  friend class fixpt;

  value_type _v;
};

// Static (arithmetic ops don't change the type). Truncation when dropping bits.
template <
  uint Sign,
  uint N_int,
  uint N_frac,
  uint Extra_flags = 0,
  class Traits     = std_fp_types_trait>
using fixpt_s
  = fixpt<Sign, N_int, N_frac, 0 | (Extra_flags & ~uint {3}), Traits>;

// Dynamic (arithmetic ops expand the type).
template <
  uint Sign,
  uint N_int,
  uint N_frac,
  uint Extra_flags = 0,
  class Traits     = std_fp_types_trait>
using fixpt_d = fixpt<
  Sign,
  N_int,
  N_frac,
  fixpt_dynamic | (Extra_flags & ~uint {3}),
  Traits>;

// Mixed (arithmetic ops expand the type until the maximum suported).
template <
  uint Sign,
  uint N_int,
  uint N_frac,
  uint Extra_flags = 0,
  class Traits     = std_fp_types_trait>
using fixpt_m
  = fixpt<Sign, N_int, N_frac, fixpt_mixed | (Extra_flags & ~uint {3}), Traits>;

//------------------------------------------------------------------------------
// convenience classes to be used as operators. The main purpose  for these is
// for float and fixed point code to be able to be geneated from the same
// source.
//------------------------------------------------------------------------------
template <int N_int, int N_frac>
struct fixpt_resize_token {};

template <
  class T,
  int N_int,
  int N_frac,
  std::enable_if_t<is_fixpt_v<T>>* = nullptr>
inline constexpr auto operator& (
  T lhs,
  fixpt_resize_token<N_int, N_frac>) noexcept
{
  return lhs.template resize<N_int, N_frac>();
}

template <int N_int1, int N_frac1, int N_int2, int N_frac2>
inline constexpr auto operator& (
  fixpt_resize_token<N_int1, N_frac1>,
  fixpt_resize_token<N_int2, N_frac2>) noexcept
{
  return fixpt_resize_token<N_int1 + N_int2, N_frac1 + N_frac2> {};
}

template <
  class T,
  int N_int,
  int N_frac,
  std::enable_if_t<!is_fixpt_v<T>>* = nullptr>
inline constexpr auto operator& (
  T lhs,
  fixpt_resize_token<N_int, N_frac>) noexcept
{
  return lhs; // passthrough for other types.
}

//------------------------------------------------------------------------------
// "ratio" related classes and functions
//------------------------------------------------------------------------------
// Create a fixed point value from a fraction.
template <
  std::intmax_t Num,
  std::intmax_t Den,
  class Traits = std_fp_types_trait>
static constexpr auto fixpt_from_ratio()
{
  constexpr auto fp       = ratio<Num, Den>::get_fixpt();
  constexpr uint max_bits = mp11::mp_back<typename Traits::conversions>::n_bits;
  constexpr uint fixed_bits = fp.n_sign - fp.n_int;
  constexpr uint max_frac   = max_bits - fixed_bits;
  constexpr uint frac_drop  = (fp.n_frac > max_frac) ? fp.n_frac - max_frac : 0;
  constexpr uint n_frac     = fp.n_frac - frac_drop;

  using fixpt_t = fixpt<fp.n_sign, fp.n_int, n_frac, fixpt_dynamic, Traits>;
  return fixpt_t::from (fp.value >> frac_drop);
}

// A class to represent a ratio truncated to some amount of bits. When used in
// arithmetic operations mixed with fixed point types it try to generate a fixed
// point constant with at most Max_frac_bits amount of fractional bits.
//
// When using "std::ratio" or "artv::ratio" constants are truncated to at most
// having the same number of fractional bits than the other operand.
template <
  std::intmax_t Num,
  std::intmax_t Den,
  uint          Max_frac_bits = (sizeof (uint) * 8)>
struct ratio_fracb {
  using base                          = ::artv::ratio<Num, Den>;
  static constexpr uint max_frac_bits = Max_frac_bits;
  static constexpr auto num           = base::num;
  static constexpr auto den           = base::den;

  template <std::intmax_t Num_v, std::intmax_t Den_v>
  using rebind = ratio_fracb<Num_v, Den_v, max_frac_bits>;
};

template <uint Max_frac_bits>
struct ratio_add_max_frac_bits_token {};

template <
  uint N_frac,
  class Ratio,
  std::enable_if_t<is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator& (Ratio, ratio_add_max_frac_bits_token<N_frac>) noexcept
{
  return ratio_fracb<Ratio::num, Ratio::den, N_frac> {};
}

// operator for adding fractional bits to a ratio, e.g:
// > (1_r/3_r) & 6_ratio_frac_bits

template <char... Chars>
constexpr auto operator"" _max_ratio_fracb()
{
  using str            = comptime::str<Chars...>;
  constexpr auto value = comptime::atoi<str>();
  static_assert (value < 1024, "_max_ratio_fracb: Bug?");
  return ratio_add_max_frac_bits_token<value> {};
}

//------------------------------------------------------------------------------
// Operators for ratio<num, den>, std::ratio<num, den>, ratio_fracb<num, den,
// max_frac_bits>
//------------------------------------------------------------------------------
namespace detail {
template <uint Max_frac_bits, std::intmax_t Num, std::intmax_t Den, class Fixpt>
static constexpr auto fixpt_from_ratio_impl (ratio<Num, Den>, Fixpt)
{
  constexpr auto fixdyn = fixpt_from_ratio<Num, Den, typename Fixpt::traits>();
  if constexpr (Fixpt::is_dynamic) {
    constexpr auto fix = fixdyn.template flags_cast<Fixpt::flags>();
    using T            = decltype (fix);
    if constexpr (T::n_frac > Max_frac_bits) {
      // limit to the same fractional size than its operand
      constexpr int n_reduce = T::n_frac - Max_frac_bits;
      return fix.template resize<0, -n_reduce>();
    }
    else {
      return fix;
    }
  }
  else {
    // The sum of bits (sign + int + frac) of static types has to match the sum
    // of bits of its underlying scalar representation. So trying to cast
    // the spec from "fixdyn" can fail on a static_assertion that verifies that.
    //
    // As the other operator is not a dynamic type we don't need to optimize
    // bits, just cat the optimized/packed ratio to the other operand's type.
    return fixdyn.template cast<Fixpt>();
  }
}

template <std::intmax_t Num, std::intmax_t Den, class Fixpt>
static constexpr uint get_max_frac_bits (::artv::ratio<Num, Den>, Fixpt)
{
  return Fixpt::n_frac;
}

template <std::intmax_t Num, std::intmax_t Den, class Fixpt>
static constexpr uint get_max_frac_bits (std::ratio<Num, Den>, Fixpt)
{
  return Fixpt::n_frac;
}

template <std::intmax_t Num, std::intmax_t Den, uint Max_frac, class Fixpt>
static constexpr uint get_max_frac_bits (ratio_fracb<Num, Den, Max_frac>, Fixpt)
{
  return Max_frac;
}

template <std::intmax_t Num, std::intmax_t Den, class Fixpt>
static constexpr auto fixpt_from_ratio (::artv::ratio<Num, Den>, Fixpt)
{
  constexpr auto r       = ::artv::ratio<Num, Den> {};
  constexpr auto f       = Fixpt {};
  constexpr auto maxfrac = get_max_frac_bits (r, f);
  return fixpt_from_ratio_impl<maxfrac> (r, f);
}

template <std::intmax_t Num, std::intmax_t Den, class Fixpt>
static constexpr auto fixpt_from_ratio (std::ratio<Num, Den>, Fixpt f)
{
  return fixpt_from_ratio (artv::ratio<Num, Den> {}, f);
}

template <std::intmax_t Num, std::intmax_t Den, uint Max_frac, class Fixpt>
static constexpr auto fixpt_from_ratio (ratio_fracb<Num, Den, Max_frac>, Fixpt)
{
  constexpr auto r        = ratio_fracb<Num, Den, Max_frac> {};
  constexpr auto f        = Fixpt {};
  constexpr uint fracbits = get_max_frac_bits (r, f);
  return fixpt_from_ratio_impl<fracbits> (artv::ratio<Num, Den> {}, f);
}

template <bool Invert, class T, class Ratio>
constexpr auto ratio_mul (T lhs, Ratio rhs) noexcept
{
  constexpr auto num = Invert ? Ratio::den : Ratio::num;
  constexpr auto den = Invert ? Ratio::num : Ratio::den;
  using rat          = ::artv::ratio<num, den>;
  constexpr auto v   = rat::get_fixpt();

  constexpr uint maxfrac = get_max_frac_bits (Ratio {}, T {});

  constexpr auto v_abs = v.n_sign ? -v.value : v.value;
  if constexpr (is_pow2 (v_abs)) {
    // for powers of two the int or fractional bits tell the shift amount
    static_assert ((v.n_int && !v.n_frac) || (!v.n_int && v.n_frac));
    // correct that multiplying by 1 causes no shift
    constexpr int int_shift = (v.n_int > 0) ? (int) v.n_int - 1 : 0;
    return lhs.template ash<int_shift - ((int) v.n_frac)>();
  }
  else if constexpr (lhs.is_dynamic) {
    // check if integer or fractional bits can be dropped.
    if constexpr (v.n_int == 0) {
      // integer bits drop (move right)
      constexpr int rshift   = v.n_frac - last_bit_set (v_abs);
      constexpr int int_drop = -std::min (rshift, (int) T::n_int);
      return (lhs * fixpt_from_ratio_impl<maxfrac> (rat {}, lhs))
        .template resize<int_drop>();
    }
    else if constexpr (v.n_frac == 0) {
      // fractional bits drop (move left)
      constexpr int lshift    = first_bit_set (v_abs) - 1;
      constexpr int frac_drop = -std::min (lshift, (int) T::n_frac);
      return (lhs * fixpt_from_ratio_impl<maxfrac> (rat {}, lhs))
        .template resize<0, frac_drop>();
    }
    else {
      return lhs * fixpt_from_ratio_impl<maxfrac> (rat {}, lhs);
    }
  }
  else {
    return lhs * fixpt_from_ratio_impl<maxfrac> (rat {}, lhs);
  }
}
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator+ (T lhs, Ratio rhs) noexcept
{
  return lhs + detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator+ (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) + rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator- (T lhs, Ratio rhs) noexcept
{
  return lhs - detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator- (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) - rhs;
}
//-----------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator* (T lhs, Ratio rhs) noexcept
{
  return detail::ratio_mul<false> (lhs, rhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator* (Ratio lhs, T rhs) noexcept
{
  return rhs * lhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator/ (T lhs, Ratio rhs) noexcept
{
  return detail::ratio_mul<true> (lhs, rhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr auto operator/ (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) / rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator== (T lhs, Ratio rhs) noexcept
{
  return lhs == detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator== (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) == rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator!= (T lhs, Ratio rhs) noexcept
{
  return lhs != detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator!= (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) != rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator> (T lhs, Ratio rhs) noexcept
{
  return lhs > detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator> (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) > rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator>= (T lhs, Ratio rhs) noexcept
{
  return lhs >= detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator>= (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) >= rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator<(T lhs, Ratio rhs) noexcept
{
  return lhs < detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator<(Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) < rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator<= (T lhs, Ratio rhs) noexcept
{
  return lhs <= detail::fixpt_from_ratio (rhs, lhs);
}

template <
  class T,
  class Ratio,
  std::enable_if_t<is_fixpt_v<T> && is_ratio_v<Ratio>>* = nullptr>
constexpr bool operator<= (Ratio lhs, T rhs) noexcept
{
  return detail::fixpt_from_ratio (lhs, rhs) <= rhs;
}

// avoid circular dependencies, as assignment overload can be a free function
template <class Ratio, std::enable_if_t<is_ratio_v<Ratio>>* = nullptr>
constexpr auto to_fixpt (Ratio) noexcept
{
  return fixpt_from_ratio<Ratio::num, Ratio::den>();
}

//------------------------------------------------------------------------------
// Operators for num<T>
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator+ (T lhs, num<U> rhs) noexcept
{
  return lhs + T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator+ (num<T> lhs, U rhs) noexcept
{
  return U {lhs} + rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator- (T lhs, num<U> rhs) noexcept
{
  return lhs - T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator- (num<T> lhs, U rhs) noexcept
{
  return U {lhs} - rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator* (T lhs, num<U> rhs) noexcept
{
  return lhs * T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator* (num<T> lhs, U rhs) noexcept
{
  return U {lhs} * rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator/ (T lhs, num<U> rhs) noexcept
{
  return lhs / T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator/ (num<T> lhs, U rhs) noexcept
{
  return U {lhs} / rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator== (T lhs, num<U> rhs) noexcept
{
  return lhs == T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator== (num<T> lhs, U rhs) noexcept
{
  return U {lhs} == rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator!= (T lhs, num<U> rhs) noexcept
{
  return lhs != T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator!= (num<T> lhs, U rhs) noexcept
{
  return U {lhs} != rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator<= (T lhs, num<U> rhs) noexcept
{
  return lhs <= T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator<= (num<T> lhs, U rhs) noexcept
{
  return U {lhs} <= rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator>= (T lhs, num<U> rhs) noexcept
{
  return lhs >= T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator>= (num<T> lhs, U rhs) noexcept
{
  return U {lhs} >= rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator<(T lhs, num<U> rhs) noexcept
{
  return lhs < T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator<(num<T> lhs, U rhs) noexcept
{
  return U {lhs} < rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator> (T lhs, num<U> rhs) noexcept
{
  return lhs > T {rhs};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator> (num<T> lhs, U rhs) noexcept
{
  return U {lhs} > rhs;
}
//------------------------------------------------------------------------------
// Operators for plain vector/arithmetic types
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator+ (T lhs, U rhs) noexcept
{
  return lhs + T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator+ (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} + rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator- (T lhs, U rhs) noexcept
{
  return lhs - T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator- (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} - rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator* (T lhs, U rhs) noexcept
{
  return lhs * T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator* (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} * rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr auto operator/ (T lhs, U rhs) noexcept
{
  return lhs / T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr auto operator/ (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} / rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator== (T lhs, U rhs) noexcept
{
  return lhs == T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator== (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} == rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator!= (T lhs, U rhs) noexcept
{
  return lhs != T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator!= (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} != rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator<= (T lhs, U rhs) noexcept
{
  return lhs <= T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator<= (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} <= rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator>= (T lhs, U rhs) noexcept
{
  return lhs >= T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator>= (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} >= rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator<(T lhs, U rhs) noexcept
{
  return lhs < T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator<(T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} < rhs;
}
//------------------------------------------------------------------------------
template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<T> && is_vec_or_scalar_v<U>>* = nullptr>
constexpr bool operator> (T lhs, U rhs) noexcept
{
  return lhs > T {rhs, detail::fixpt_arith_tag {}};
}

template <
  class T,
  class U,
  std::enable_if_t<is_fixpt_v<U> && is_vec_or_scalar_v<T>>* = nullptr>
constexpr bool operator> (T lhs, U rhs) noexcept
{
  return U {lhs, detail::fixpt_arith_tag {}} > rhs;
}
//------------------------------------------------------------------------------
// ternary operators can't be overloaded. This is mostly done for native
// vector types.
template <class C, class U, std::enable_if_t<is_fixpt_v<U>>* = nullptr>
constexpr auto fixpt_select (C cond, U a, U b) noexcept
{
  return U::from (cond ? a.value() : b.value());
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
constexpr auto fixpt_abs (T v) noexcept
{
  if constexpr (v.is_signed) {
    constexpr uint shift = v.scalar_bits - 1;

    auto raw  = v.value();
    auto mask = raw >> shift;
    return T::from ((mask + raw) ^ mask);
  }
  else {
    return v;
  }
}
//------------------------------------------------------------------------------
template <class T, class U, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
constexpr auto fixpt_clamp (T v, U min, U max) noexcept
{
  assert (min <= max);
  auto kmin = (T) min;
  auto kmax = (T) max;
  auto ret  = fixpt_select (v < kmin, kmin, v);
  ret       = fixpt_select (v > kmax, kmax, v);
  return ret;
}
//------------------------------------------------------------------------------

// TODO: fixpt_ceil, fixpt_floor, fixpt_round, fixpt_exp2...
//------------------------------------------------------------------------------
} // namespace artv
