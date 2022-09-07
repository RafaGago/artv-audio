#pragma once

//------------------------------------------------------------------------------
// All the classes under parts are coded stateless for $REASONS (see README).
// These are wrappers to make most of them classes.
//------------------------------------------------------------------------------

#include <array>
#include <cassert>
#include <utility>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

namespace detail {
//------------------------------------------------------------------------------
// reminder: co = coefficients, co_int = coefficients internal, st = states.
template <
  class T_co,
  uint N_co,
  uint N_co_int,
  uint N_co_layers,
  class T_st,
  uint N_st,
  uint N_st_layers,
  bool dynamic>
class part_memory;
//------------------------------------------------------------------------------
// array-backed specialization
template <
  class T_co,
  uint N_co,
  uint N_co_int,
  uint N_co_layers,
  class T_st,
  uint N_st,
  uint N_st_layers>
class part_memory<
  T_co,
  N_co,
  N_co_int,
  N_co_layers,
  T_st,
  N_st,
  N_st_layers,
  false> {
public:
  //----------------------------------------------------------------------------
  void reset() { memset (this, 0, sizeof *this); }
  //----------------------------------------------------------------------------
  crange<T_co> get_coeffs (uint layer = 0)
  {
    assert (layer < N_co_layers);
    return _coeffs[layer];
  }
  //----------------------------------------------------------------------------
  crange<T_co> get_all_coeffs()
  {
    return {_coeffs[0].data(), _coeffs.size() * _coeffs[0].size()};
  }
  //----------------------------------------------------------------------------
  crange<T_co> get_coeffs_int (uint layer = 0)
  {
    assert (layer < N_co_layers);
    return _coeffs_int[layer];
  }
  //----------------------------------------------------------------------------
  crange<T_co> get_all_coeffs_int()
  {
    return {_coeffs_int[0].data(), _coeffs_int.size() * _coeffs_int[0].size()};
  }
  //----------------------------------------------------------------------------
  crange<T_st> get_states (uint layer = 0)
  {
    assert (layer < N_st_layers);
    return _states[layer];
  }
  //----------------------------------------------------------------------------
  crange<T_st> get_all_states()
  {
    return {_states[0].data(), _states.size() * _states[0].size()};
  }
  //----------------------------------------------------------------------------
private:
  alignas (T_co) std::array<std::array<T_co, N_co>, N_co_layers> _coeffs;
  alignas (
    T_co) std::array<std::array<T_co, N_co_int>, N_co_layers> _coeffs_int;
  alignas (T_st) std::array<std::array<T_st, N_st>, N_st_layers> _states;
};
//------------------------------------------------------------------------------
// vector-backed specialization, same types
// reminder: co = coefficients, co_int = coefficients internal, st = states.
template <
  class T,
  uint N_co,
  uint N_co_int,
  uint N_co_layers,
  uint N_st,
  uint N_st_layers>
class part_memory<T, N_co, N_co_int, N_co_layers, T, N_st, N_st_layers, true> {
public:
  //----------------------------------------------------------------------------
  void reset()
  {
    _mem.clear();
    _mem.resize ((N_co + N_co_int) * N_co_layers + (N_st * N_st_layers));
  }
  //----------------------------------------------------------------------------
  crange<T> get_coeffs (uint layer = 0)
  {
    assert (layer < N_co_layers);
    return {&_mem[coeff_offset + layer * N_co], N_co};
  }
  //----------------------------------------------------------------------------
  crange<T> get_all_coeffs()
  {
    return {&_mem[coeff_offset], N_co * N_co_layers};
  }
  //----------------------------------------------------------------------------
  crange<T> get_coeffs_int (uint layer = 0)
  {
    assert (layer < N_co_layers);
    return {&_mem[coeff_int_offset + layer * N_co_int], N_co_int};
  }
  //----------------------------------------------------------------------------
  crange<T> get_all_coeffs_int()
  {
    return {&_mem[coeff_int_offset], N_co_int * N_co_layers};
  }
  //----------------------------------------------------------------------------
  crange<T> get_states (uint layer = 0)
  {
    assert (layer < N_st_layers);
    return {&_mem[states_offset + layer * N_st], N_st};
  }
  //----------------------------------------------------------------------------
  crange<T> get_all_states()
  {
    return {&_mem[states_offset], N_st * N_st_layers};
  }
  //----------------------------------------------------------------------------
private:
  static constexpr uint coeff_offset     = 0;
  static constexpr uint coeff_int_offset = N_co * N_co_layers;
  static constexpr uint states_offset
    = coeff_int_offset + (N_co_int * N_co_layers);

  std::vector<T, overaligned_allocator<T, sizeof (T)>> _mem;
};
//------------------------------------------------------------------------------
// vector-backed specialization, different types, could be the only one used
// instead of the same-type variant, which was done for the debugger's
// convenience.
// reminder: co = coefficients, co_int = coefficients internal, st = states.
template <
  class T_co,
  uint N_co,
  uint N_co_int,
  uint N_co_layers,
  class T_st,
  uint N_st,
  uint N_st_layers>
class part_memory<
  T_co,
  N_co,
  N_co_int,
  N_co_layers,
  T_st,
  N_st,
  N_st_layers,
  true> {
public:
  //----------------------------------------------------------------------------
  void reset()
  {
    _mem.clear();
    _mem.resize (n_bytes);
  }
  //----------------------------------------------------------------------------
  crange<T_co> get_coeffs (uint layer = 0)
  {
    assert (layer < N_co_layers);
    // This is intended to be used with __may_alias__ types, hence the C casts
    auto offset = (T_co*) &_mem[co_offset];
    return {offset + (layer * N_co), N_co};
  }
  //----------------------------------------------------------------------------
  crange<T_co> get_all_coeffs()
  {
    // This is intended to be used with __may_alias__ types, hence the C casts
    auto offset = (T_co*) &_mem[co_offset];
    return {offset, N_co * N_co_layers};
  }
  //----------------------------------------------------------------------------
  crange<T_co> get_coeffs_int (uint layer = 0)
  {
    assert (layer < N_co_layers);
    // This is intended to be used with __may_alias__ types, hence the C casts
    auto offset = (T_co*) &_mem[co_int_offset];
    return {offset + (layer * N_co_int), N_co_int};
  }
  //----------------------------------------------------------------------------
  crange<T_co> get_all_coeffs_int()
  {
    // This is intended to be used with __may_alias__ types, hence the C casts
    auto offset = (T_co*) &_mem[co_int_offset];
    return {offset, N_co_int * N_co_layers};
  }
  //----------------------------------------------------------------------------
  crange<T_st> get_states (uint layer = 0)
  {
    // This is intended to be used with __may_alias__ types, hence the C casts
    auto offset = (T_st*) &_mem[st_offset];
    return {offset + (layer * N_st), N_st};
  }
  //----------------------------------------------------------------------------
  crange<T_st> get_all_states()
  {
    // This is intended to be used with __may_alias__ types, hence the C casts
    auto offset = (T_st*) &_mem[st_offset];
    return {offset, N_st * N_st_layers};
  }
  //----------------------------------------------------------------------------
private:
  // this could work generically by looking at the first bit set of the size of
  // the type, but it is not portable and constexpr until C++20
  // (std::countr_zero). Restricting to arithmetic types, which are powers of
  // two, so if the first is the biggest the second is guaranteed to be
  // correctly aligned.
  static_assert (std::is_arithmetic_v<T_co> || is_vec_v<T_co>, "");
  static_assert (std::is_arithmetic_v<T_st> || is_vec_v<T_co>, "");

  static constexpr bool co_is_big_type = sizeof (T_co) >= sizeof (T_st);
  using big_type   = std::conditional_t<co_is_big_type, T_co, T_st>;
  using small_type = std::conditional_t<!co_is_big_type, T_co, T_st>;

  static constexpr uint n_co_bytes     = N_co * N_co_layers * sizeof (T_co);
  static constexpr uint n_co_int_bytes = N_co_int * N_co_layers * sizeof (T_co);
  static constexpr uint n_st_bytes     = N_st * N_st_layers * sizeof (T_st);
  static constexpr uint n_bytes = n_co_bytes + n_co_int_bytes + n_st_bytes;

  static constexpr uint co_offset     = co_is_big_type ? 0 : n_st_bytes;
  static constexpr uint co_int_offset = co_offset + n_co_bytes;
  static constexpr uint st_offset
    = co_is_big_type ? (co_int_offset + n_co_int_bytes) : 0;

  std::vector<u8, overaligned_allocator<u8, sizeof (big_type)>> _mem;
};
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
// A DSP part to hold an array of instances of the same type.
template <class Part, class Vect, uint Size = 1, bool dynamic_mem = false>
struct part_class_array {
public:
  static_assert (is_vec_of_float_type_v<Vect>, "");
  static_assert (Size > 0, "");
  //----------------------------------------------------------------------------
  using part                     = Part;
  using value_type               = Vect;
  using builtin                  = vec_value_type_t<value_type>;
  static constexpr uint size     = Size;
  static constexpr uint n_coeffs = part::n_coeffs;
  //----------------------------------------------------------------------------
  void reset() { _mem.reset(); }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <class... Ts>
  void reset_coeffs_ext (uint idx, crange<value_type> coeff_out, Ts&&... args)
  {
    assert (coeff_out.size() >= part::n_coeffs);

    if constexpr (part::n_coeffs_int == 0) {
      part::template reset_coeffs<value_type> (
        coeff_out, std::forward<Ts> (args)...);
    }
    else {
      auto co_int = _mem.get_coeffs_int (idx);
      part::template reset_coeffs<value_type> (
        coeff_out, co_int, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs_on_idx (uint idx, Ts&&... args)
  {
    auto co = get_coeffs (idx);
    reset_coeffs_ext (idx, co, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_states_on_idx (uint idx, Ts&&... args)
  {
    auto st = _mem.get_states (idx);
    part::template reset_states<value_type> (st, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    auto co = get_coeffs (idx);
    auto st = _mem.get_states (idx);

    if constexpr (part::n_coeffs_int == 0) {
      return part::template tick<value_type> (
        co, st, std::forward<Ts> (args)...);
    }
    else {
      auto co_int = _mem.get_coeffs_int (idx);
      return part::template tick<value_type> (
        co, co_int, st, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    reset_coeffs_on_idx (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_states (Ts&&... args)
  {
    reset_states_on_idx (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick (Ts&&... args)
  {
    return tick_on_idx (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    static_assert (Idx < size);
    reset_coeffs_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void reset_states (Ts&&... args)
  {
    static_assert (Idx < size);
    reset_states_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < size);
    return tick_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs_cascade (Ts&&... args)
  {
    // This might not compile if all internal FX don't use the same paramters
    for (uint i = 0; i < size; ++i) {
      reset_coeffs_on_idx (i, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_states_cascade (Ts&&... args)
  {
    for (uint i = 0; i < size; ++i) {
      reset_states_on_idx (i, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class... Ts>
  auto tick_cascade (T in, Ts&&... args)
  {
    // This might not compile if all internal FX doesn't use the same parameters
    // and return types
    T out = in;
    for (uint i = 0; i < size; ++i) {
      out = tick_on_idx (i, out, std::forward<Ts> (args)...);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<value_type> get_coeffs (uint idx = 0) { return _mem.get_coeffs (idx); }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<value_type> get_all_coeffs() { return _mem.get_all_coeffs(); }
  //----------------------------------------------------------------------------
  // for initial conditions setup
  crange<value_type> get_states (uint idx = 0) { return _mem.get_states (idx); }
  //----------------------------------------------------------------------------
  void zero_all_states() { crange_memset (_mem.get_all_states(), 0); }
  //----------------------------------------------------------------------------
private:
  detail::part_memory<
    value_type,
    part::n_coeffs,
    part::n_coeffs_int,
    size,
    value_type,
    part::n_states,
    size,
    dynamic_mem>
    _mem;
};
//------------------------------------------------------------------------------
// A DSP part to (optionally) an array of instances that uses a single
// coefficient set for element of the array and is of width = 1 independently of
// the passed vector (Vect) width. Useful for things like e.g. DC blockers
template <class Part, class Vect, uint Size = 1, bool dynamic_mem = false>
struct part_class_array_coeffs_global {
public:
  static_assert (is_vec_of_float_type_v<Vect>, "");
  static_assert (Size > 0, "");
  //----------------------------------------------------------------------------
  using part                     = Part;
  using value_type               = Vect;
  using builtin                  = vec_value_type_t<value_type>;
  using co_vec_type              = vec<builtin, 1>;
  static constexpr uint n_coeffs = part::n_coeffs;
  static constexpr uint size     = Size;
  //----------------------------------------------------------------------------
  void reset() { _mem.reset(); }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <class... Ts>
  void reset_coeffs_ext (crange<builtin> coeff_out, Ts&&... args)
  {
    assert (coeff_out.size() >= part::n_coeffs);

    if constexpr (part::n_coeffs_int == 0) {
      part::template reset_coeffs<co_vec_type> (
        coeff_out.cast (co_vec_type {}), std::forward<Ts> (args)...);
    }
    else {
      auto co_int = _mem.get_coeffs_int();
      part::template reset_coeffs<co_vec_type> (
        coeff_out.cast (co_vec_type {}),
        co_int.cast (co_vec_type {}),
        std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    auto co = get_coeffs();
    reset_coeffs_ext (co, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_states_on_idx (uint idx, Ts&&... args)
  {
    auto st = _mem.get_states (idx);
    part::template reset_states<value_type> (st, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    auto co = get_coeffs();
    auto st = _mem.get_states (idx);

    if constexpr (part::n_coeffs_int == 0) {
      return part::template tick<value_type> (
        co, st, std::forward<Ts> (args)...);
    }
    else {
      auto co_int = _mem.get_coeffs_int();
      return part::template tick<value_type> (
        co, co_int, st, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_states (Ts&&... args)
  {
    reset_states_on_idx (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick (Ts&&... args)
  {
    return tick_on_idx (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void reset_states (Ts&&... args)
  {
    static_assert (Idx < size);
    reset_states_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < size);
    return tick_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_states_cascade (Ts&&... args)
  {
    for (uint i = 0; i < size; ++i) {
      reset_states_on_idx (i, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class... Ts>
  auto tick_cascade (T in, Ts&&... args)
  {
    // This might not compile if all internal FX doesn't use the same parameters
    // and return types
    T out = in;
    for (uint i = 0; i < size; ++i) {
      out = tick_on_idx (i, out, std::forward<Ts> (args)...);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  crange<builtin> get_coeffs() { return _mem.get_coeffs(); }
  crange<builtin> get_all_coeffs() { return _mem.get_all_coeffs(); }
  //----------------------------------------------------------------------------
  crange<value_type> get_states (uint idx = 0) { return _mem.get_states (idx); }
  void zero_all_states() { crange_memset (_mem.get_all_states(), 0); }
  //----------------------------------------------------------------------------
private:
  detail::part_memory<
    builtin,
    part::n_coeffs,
    part::n_coeffs_int,
    1,
    value_type,
    part::n_states,
    size,
    dynamic_mem>
    _mem;
};
//------------------------------------------------------------------------------
template <class PartList, class Vect, bool dynamic_mem = false>
struct part_classes;

// Different DSP parts stored on a single class and accessed by compile-time
// index.
template <
  template <class...>
  class List,
  class Vect,
  class... Parts,
  bool dynamic_mem>
struct part_classes<List<Parts...>, Vect, dynamic_mem> {
public:
  static_assert (sizeof...(Parts) > 0);
  static_assert (is_vec_of_float_type_v<Vect>, "");
  //----------------------------------------------------------------------------
  using parts                   = mp_list<Parts...>;
  using value_type              = Vect;
  using builtin                 = vec_value_type_t<value_type>;
  static constexpr uint n_parts = sizeof...(Parts);

  template <uint Idx>
  using get_part = mp11::mp_at_c<parts, Idx>;
  //----------------------------------------------------------------------------
  void reset() { _mem.reset(); }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <uint Idx, class... Ts>
  void reset_coeffs_ext (crange<value_type> coeff_out, Ts&&... args)
  {
    static_assert (Idx < sizeof...(Parts));
    using part = get_part<Idx>;

    assert (coeff_out.size() >= part::n_coeffs);

    if constexpr (part::n_coeffs_int == 0) {
      part::template reset_coeffs<value_type> (
        coeff_out, std::forward<Ts> (args)...);
    }
    else {
      auto co_int = get_coeffs_int<Idx, part::n_coeffs_int>();
      part::template reset_coeffs<value_type> (
        coeff_out, co_int, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    static_assert (Idx < sizeof...(Parts));
    using part = get_part<Idx>;

    auto co = get_coeffs<Idx, part::n_coeffs>();
    part::template reset_coeffs<value_type> (co, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void reset_states (Ts&&... args)
  {
    static_assert (Idx < sizeof...(Parts));
    using part = get_part<Idx>;

    auto st = get_states<Idx, part::n_states>();
    part::template reset_states<value_type> (st, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < sizeof...(Parts));
    using part = get_part<Idx>;

    auto co = get_coeffs<Idx, part::n_coeffs>();
    auto st = get_states<Idx, part::n_states>();

    if constexpr (part::n_coeffs_int == 0) {
      return part::template tick<value_type> (
        co, st, std::forward<Ts> (args)...);
    }
    else {
      auto co_int = get_coeffs_int<Idx, part::n_coeffs_int>();
      return part::template tick<value_type> (
        co, co_int, st, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick (Ts&&... args)
  {
    return tick<0> (std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs_cascade (Ts&&... args)
  {
    mp11::mp_for_each<mp11::mp_iota_c<sizeof...(Parts)>> ([&, this] (auto val) {
      static constexpr uint idx = decltype (val)::value;
      reset_coeffs<idx> (std::forward<Ts> (args)...);
    });
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_states_cascade (Ts&&... args)
  {
    mp11::mp_for_each<mp11::mp_iota_c<sizeof...(Parts)>> ([&, this] (auto val) {
      static constexpr uint idx = decltype (val)::value;
      reset_states<idx> (std::forward<Ts> (args)...);
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class... Ts>
  T tick_cascade (T in, Ts&&... args)
  {
    // This might not compile if all internal FX doesn't use the same
    // parameters and return types
    T out = in;
    mp11::mp_for_each<mp11::mp_iota_c<sizeof...(Parts)>> ([&, this] (auto val) {
      static constexpr uint idx = decltype (val)::value;
      out                       = tick<idx> (out, std::forward<Ts> (args)...);
    });
    return out;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  crange<value_type> get_coeffs()
  {
    static_assert (Idx < sizeof...(Parts));
    using part = get_part<Idx>;
    return get_coeffs<Idx, part::n_coeffs>();
  }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<value_type> get_all_coeffs() { return _mem.get_all_coeffs(); }
  //----------------------------------------------------------------------------
  template <uint Idx>
  static constexpr uint get_coeff_offset()
  {
    static_assert (Idx < n_parts, "");
    return get_offset<n_coeffs_list, Idx>();
  }
  //----------------------------------------------------------------------------
  void zero_all_states() { crange_memset (_mem.get_all_states(), 0); }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <uint Idx>
  crange<value_type> get_states()
  {
    using part = get_part<Idx>;
    return get_states<Idx, part::n_states>();
  }
  //----------------------------------------------------------------------------
  template <uint Idx, uint N_states>
  crange<value_type> get_states()
  {
    if constexpr (N_states > 0) {
      constexpr auto offset = get_offset<n_states_list, Idx>();
      auto           st     = _mem.get_states (0);
      return {&st[offset], N_states};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, uint N_coeffs>
  crange<value_type> get_coeffs_int()
  {
    if constexpr (N_coeffs > 0) {
      constexpr auto offset = get_offset<n_coeffs_int_list, Idx>();

      auto co = _mem.get_coeffs_int (0);
      return {&co[offset], N_coeffs};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, uint N_coeffs>
  crange<value_type> get_coeffs()
  {
    if constexpr (N_coeffs > 0) {
      constexpr auto offset = get_coeff_offset<Idx>();

      auto co = _mem.get_coeffs (0);
      return {&co[offset], N_coeffs};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  template <class SizeList, uint End>
  static constexpr uint get_offset()
  {
    // for some reason I didn't get working C++17's fold expressions with
    // compile time types. TODO (...or not, live is too short).
    using lst = mp11::mp_take_c<SizeList, End>;
    return mp11::mp_fold<lst, k_int<0>, mp11::mp_plus>::value;
  }
  //----------------------------------------------------------------------------
  template <class T>
  using to_n_coeffs = k_int<T::n_coeffs>;
  template <class T>
  using to_n_coeffs_int = k_int<T::n_coeffs_int>;
  template <class T>
  using to_n_states = k_int<T::n_states>;

  using n_coeffs_list     = mp11::mp_transform<to_n_coeffs, parts>;
  using n_coeffs_int_list = mp11::mp_transform<to_n_coeffs_int, parts>;
  using n_states_list     = mp11::mp_transform<to_n_states, parts>;

  static constexpr uint n_coeffs_total = get_offset<n_coeffs_list, n_parts>();
  static constexpr uint n_coeffs_int_total
    = get_offset<n_coeffs_int_list, n_parts>();
  static constexpr uint n_states_total = get_offset<n_states_list, n_parts>();

  detail::part_memory<
    value_type,
    n_coeffs_total,
    n_coeffs_int_total,
    1,
    value_type,
    n_states_total,
    1,
    dynamic_mem>
    _mem;

public:
  static constexpr uint n_coeffs = n_coeffs_total;
};
//------------------------------------------------------------------------------
namespace detail {

template <
  class PartList,
  class Vect,
  class CoeffType,
  uint Size,
  bool dynamic_mem>
struct parts_union_array;

template <
  template <class...>
  class List,
  class Vect,
  class CoeffType,
  uint Size,
  bool dynamic_mem,
  class... Parts>
struct parts_union_array<List<Parts...>, Vect, CoeffType, Size, dynamic_mem> {
public:
  static_assert (is_vec_of_float_type_v<Vect>, "");
  static_assert (Size > 0, "");
  //----------------------------------------------------------------------------
  using parts                   = mp_list<Parts...>;
  static constexpr uint n_parts = sizeof...(Parts);
  static_assert (n_parts > 0);
  using value_type           = Vect;
  using builtin              = vec_value_type_t<value_type>;
  using coeff_type           = CoeffType;
  static constexpr uint size = Size;

  static_assert (
    std::is_same_v<Vect, CoeffType> || std::is_same_v<builtin, CoeffType>,
    "");
  //----------------------------------------------------------------------------
  void reset() { _mem.reset(); }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <class Part, class... Ts>
  void reset_coeffs_ext (uint idx, crange<coeff_type> coeff_out, Ts&&... args)
  {
    static_assert (mp11::mp_find<parts, Part>::value < n_parts, "");

    assert (coeff_out.size() >= Part::n_coeffs);

    if constexpr (std::is_same_v<value_type, coeff_type>) {
      if constexpr (Part::n_coeffs_int == 0) {
        Part::template reset_coeffs<value_type> (
          coeff_out, std::forward<Ts> (args)...);
      }
      else {
        auto co_int = _mem.get_coeffs_int (idx);
        Part::template reset_coeffs<value_type> (
          coeff_out, co_int, std::forward<Ts> (args)...);
      }
    }
    else {
      using x1_t = vec<builtin, 1>;
      if constexpr (Part::n_coeffs_int == 0) {
        Part::template reset_coeffs<x1_t> (
          make_crange (coeff_out).cast (x1_t {}), std::forward<Ts> (args)...);
      }
      else {
        auto co_int = _mem.get_coeffs_int (idx);
        Part::template reset_coeffs<x1_t> (
          make_crange (coeff_out).cast (x1_t {}),
          co_int.cast (x1_t {}),
          std::forward<Ts> (args)...);
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  void reset_coeffs_on_idx (uint idx, Ts&&... args)
  {
    auto co = get_coeffs (idx);
    reset_coeffs_ext<Part> (idx, co, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  void reset_states_on_idx (uint idx, Ts&&... args)
  {
    static_assert (mp11::mp_find<parts, Part>::value < n_parts, "");

    auto st = _mem.get_states (idx);
    Part::template reset_states<value_type> (st, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  // generic 0 setting, for safety use the variant above
  void reset_states_on_idx (uint idx)
  {
    auto st = _mem.get_states (idx);
    crange_memset (st, 0);
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    static_assert (mp11::mp_find<parts, Part>::value < n_parts, "");

    auto co = get_coeffs (idx);
    auto st = _mem.get_states (idx);

    if constexpr (Part::n_coeffs_int == 0) {
      return Part::template tick<value_type> (
        co, st, std::forward<Ts> (args)...);
    }
    else {
      auto co_int = _mem.get_coeffs_int (idx);
      return Part::template tick<value_type> (
        co, co_int, st, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    reset_coeffs_on_idx<Part> (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  void reset_states (Ts&&... args)
  {
    reset_states_on_idx<Part> (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  auto tick (Ts&&... args)
  {
    return tick_on_idx<Part> (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, uint Idx, class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    static_assert (Idx < size);
    reset_coeffs_on_idx<Part> (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, uint Idx, class... Ts>
  void reset_states (Ts&&... args)
  {
    static_assert (Idx < size);
    reset_states_on_idx<Part> (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < size);
    return tick_on_idx<Part> (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<coeff_type> get_coeffs (uint idx = 0) { return _mem.get_coeffs (idx); }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<coeff_type> get_all_coeffs() { return _mem.get_all_coeffs(); }
  //----------------------------------------------------------------------------
  crange<coeff_type> get_states (uint idx = 0) { return _mem.get_states (idx); }
  void zero_all_states() { crange_memset (_mem.get_all_states(), 0); }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  using to_n_coeffs = k_int<T::n_coeffs>;
  template <class T>
  using to_n_coeffs_int = k_int<T::n_coeffs_int>;
  template <class T>
  using to_n_states = k_int<T::n_states>;

  using n_coeffs_types     = mp11::mp_transform<to_n_coeffs, parts>;
  using n_coeffs_int_types = mp11::mp_transform<to_n_coeffs_int, parts>;
  using n_states_types     = mp11::mp_transform<to_n_states, parts>;

  static constexpr uint max_coeffs
    = mp11::mp_max_element<n_coeffs_types, mp11::mp_less>::value;

  static constexpr uint max_coeffs_int
    = mp11::mp_max_element<n_coeffs_int_types, mp11::mp_less>::value;

  static constexpr uint max_states
    = mp11::mp_max_element<n_states_types, mp11::mp_less>::value;

  detail::part_memory<
    coeff_type,
    max_coeffs,
    max_coeffs_int,
    size,
    value_type,
    max_states,
    size,
    dynamic_mem>
    _mem;

public:
  static constexpr uint n_coeffs = max_coeffs;
};
} // namespace detail

//------------------------------------------------------------------------------
// N instances of one of many different DSP parts. Every instance can be only
// one of the available DSP parts simultaneously. This is controlled
// externally.
template <class PartList, class Vect, uint Size = 1, bool dynamic_mem = false>
using parts_union_array
  = detail::parts_union_array<PartList, Vect, Vect, Size, dynamic_mem>;
//------------------------------------------------------------------------------
// As "parts_union_array" but the coefficient set is of width = 1,
// independently of the with of "Vect"
template <class PartList, class Vect, uint Size = 1, bool dynamic_mem = false>
using parts_union_array_coeffs_by_idx = detail::
  parts_union_array<PartList, Vect, vec_value_type_t<Vect>, Size, dynamic_mem>;
//------------------------------------------------------------------------------
} // namespace artv
