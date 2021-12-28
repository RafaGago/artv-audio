#pragma once

//------------------------------------------------------------------------------
// All the classes under parts are coded stateless for $REASONS (see README).
// These are wrappers to make most of them classes.
//------------------------------------------------------------------------------

#include <array>
#include <cassert>
#include <tuple>
#include <utility>

#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// A DSP part to (optionally) an array of instances.
template <class Vect, class Part, uint Count = 1>
struct part_to_class {
public:
  static_assert (is_vec_of_float_type_v<Vect>, "");
  static_assert (Count > 0, "");
  //----------------------------------------------------------------------------
  using part                    = Part;
  using value_type              = Vect;
  using builtin                 = vec_value_type_t<value_type>;
  static constexpr uint n_elems = Count;
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <class... Ts>
  static void get_target_coeffs (crange<value_type> coeffs, Ts&&... args)
  {
    assert (coeffs.size() >= part::n_coeffs);
    part::template reset_coeffs<value_type> (
      coeffs, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs_on_idx (uint idx, Ts&&... args)
  {
    assert (idx < n_elems);
    get_target_coeffs (_coeffs[idx], std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  void reset_states_on_idx (uint idx)
  {
    assert (idx < n_elems);
    part::template reset_states<value_type> (_states[idx]);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    assert (idx < n_elems);
    return part::template tick<value_type> (
      _coeffs[idx], _states[idx], std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    reset_coeffs_on_idx (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  void reset_states() { reset_states_on_idx (0); }
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
    static_assert (Idx < n_elems);
    reset_coeffs_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void reset_states()
  {
    static_assert (Idx < n_elems);
    reset_states_on_idx (Idx);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < n_elems);
    return tick_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs_cascade (Ts&&... args)
  {
    // This might not compile if all internal FX don't use the same paramters
    for (uint i = 0; i < n_elems; ++i) {
      reset_coeffs_on_idx (i, std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  void reset_states_cascade()
  {
    for (uint i = 0; i < n_elems; ++i) {
      reset_states_on_idx (i);
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class... Ts>
  auto tick_cascade (T in, Ts&&... args)
  {
    // This might not compile if all internal FX doesn't use the same parameters
    // and return types
    T out = in;
    for (uint i = 0; i < n_elems; ++i) {
      out = tick_on_idx (i, out, std::forward<Ts> (args)...);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  std::array<value_type, part::n_coeffs>& get_coeffs()
  {
    static_assert (Idx < n_elems, "");
    return _coeffs[Idx];
  }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<value_type> get_coeffs (uint idx = 0)
  {
    assert (idx < n_elems);
    return (idx < n_elems) ? make_crange (_coeffs[idx]) : crange<value_type> {};
  }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<value_type> get_all_coeffs()
  {
    return {_coeffs.data(), _coeffs.size() * _coeffs[0].size()};
  }
  //----------------------------------------------------------------------------
private:
  alignas (sizeof (value_type))
    std::array<std::array<value_type, part::n_coeffs>, n_elems> _coeffs;
  alignas (sizeof (value_type))
    std::array<std::array<value_type, part::n_states>, n_elems> _states;
};
//------------------------------------------------------------------------------
// A DSP part to (optionally) an array of instances that uses a single
// coefficient set for element of the array and is of width = 1 independently of
// the passed vector (Vect) width. Useful for things like e.g. DC blockers
template <class Vect, class Part, uint Count = 1>
struct part_to_class_single_coeff_all {
public:
  static_assert (is_vec_of_float_type_v<Vect>, "");
  static_assert (Count > 0, "");
  //----------------------------------------------------------------------------
  using part       = Part;
  using value_type = Vect;
  using builtin    = vec_value_type_t<value_type>;

  static constexpr uint n_elems = Count;
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <class... Ts>
  static void get_target_coeffs (crange<builtin> coeffs, Ts&&... args)
  {
    assert (coeffs.size() >= part::n_coeffs);
    using x1_t = vec<builtin, 1>;
    part::template reset_coeffs<x1_t> (
      coeffs.cast (x1_t {}), std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    get_target_coeffs (_coeffs, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  void reset_states_on_idx (uint idx)
  {
    assert (idx < n_elems);
    part::template reset_states<value_type> (_states[idx]);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    assert (idx < n_elems);
    return part::template tick<value_type> (
      make_crange (_coeffs),
      make_crange (_states[idx]),
      std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  void reset_states() { reset_states_on_idx (0); }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick (Ts&&... args)
  {
    return tick_on_idx (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void reset_states()
  {
    static_assert (Idx < n_elems);
    reset_states_on_idx (Idx);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < n_elems);
    return tick_on_idx (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  void reset_states_cascade()
  {
    for (uint i = 0; i < n_elems; ++i) {
      reset_states_on_idx (i);
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class... Ts>
  auto tick_cascade (T in, Ts&&... args)
  {
    // This might not compile if all internal FX doesn't use the same parameters
    // and return types
    T out = in;
    for (uint i = 0; i < n_elems; ++i) {
      out = tick_on_idx (i, out, std::forward<Ts> (args)...);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  std::array<value_type, part::n_coeffs>& get_coeffs() { return _coeffs; }
  //----------------------------------------------------------------------------
private:
  alignas (sizeof (value_type)) std::array<builtin, part::n_coeffs> _coeffs;
  alignas (sizeof (value_type))
    std::array<std::array<value_type, part::n_states>, n_elems> _states;
};
//------------------------------------------------------------------------------
// Different DSP parts stored on a single class and accessed by index.
template <class Vect, class... Parts>
struct parts_to_class {
public:
  static_assert (sizeof...(Parts) > 0);
  static_assert (is_vec_of_float_type_v<Vect>, "");
  //----------------------------------------------------------------------------
  using parts      = mp_list<Parts...>;
  using value_type = Vect;
  using builtin    = vec_value_type_t<value_type>;
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <uint Idx, class... Ts>
  static void get_target_coeffs (crange<value_type> coeffs, Ts&&... args)
  {
    static_assert (Idx < sizeof...(Parts));
    using part = get_part<Idx>;
    assert (coeffs.size() >= part::n_coeffs);

    part::template reset_coeffs<value_type> (
      coeffs, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    auto& coeffs = std::get<Idx> (_coeffs);
    get_target_coeffs<Idx> (coeffs, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void reset_states()
  {
    static_assert (Idx < sizeof...(Parts));
    using part   = get_part<Idx>;
    auto& states = std::get<Idx> (_states);
    part::template reset_states<value_type> (_states);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < sizeof...(Parts));
    using part   = get_part<Idx>;
    auto& coeffs = std::get<Idx> (_coeffs);
    auto& states = std::get<Idx> (_states);
    return part::template tick<value_type> (
      _coeffs, _states, std::forward<Ts> (args)...);
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
  void reset_states_cascade()
  {
    mp11::mp_for_each<mp11::mp_iota_c<sizeof...(Parts)>> ([&, this] (auto val) {
      static constexpr uint idx = decltype (val)::value;
      reset_states<idx>();
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class... Ts>
  auto tick_cascade (T in, Ts&&... args)
  {
    // This might not compile if all internal FX doesn't use the same parameters
    // and return types
    T out = in;
    mp11::mp_for_each<mp11::mp_iota_c<sizeof...(Parts)>> ([&, this] (auto val) {
      static constexpr uint idx = decltype (val)::value;
      out                       = tick<idx> (out, std::forward<Ts> (args)...);
    });
    return out;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  auto& get_coeffs()
  {
    static_assert (Idx < sizeof...(Parts));
    return std::get<Idx> (_coeffs);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  using coeffs_array = std::array<Vect, T::n_coeffs>;
  template <class T>
  using states_array = std::array<Vect, T::n_states>;

  template <uint Idx>
  using get_part = mp11::mp_at_c<mp_list<Parts...>, Idx>;

  alignas (sizeof (value_type)) std::tuple<coeffs_array<Parts>...> _coeffs;
  alignas (sizeof (value_type)) std::tuple<states_array<Parts>...> _states;
};
//------------------------------------------------------------------------------
namespace detail {
template <class Vect, class CoeffType, uint Count, class... Parts>
struct parts_to_class_one_of {
public:
  static_assert (is_vec_of_float_type_v<Vect>, "");
  static_assert (Count > 0, "");
  //----------------------------------------------------------------------------
  using parts                   = mp_list<Parts...>;
  static constexpr uint n_parts = sizeof...(Parts);
  static_assert (n_parts > 0);
  using value_type              = Vect;
  using builtin                 = vec_value_type_t<value_type>;
  using coeff_type              = CoeffType;
  static constexpr uint n_elems = Count;

  static_assert (
    std::is_same_v<Vect, CoeffType> || std::is_same_v<builtin, CoeffType>,
    "");
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  template <class Part, class... Ts>
  static void get_target_coeffs (crange<coeff_type> coeffs, Ts&&... args)
  {
    static_assert (mp11::mp_find<parts, Part>::value < n_parts, "");

    if constexpr (std::is_same_v<value_type, coeff_type>) {
      Part::template reset_coeffs<value_type> (
        coeffs, std::forward<Ts> (args)...);
    }
    else {
      using x1_t = vec<builtin, 1>;
      Part::template reset_coeffs<x1_t> (
        make_crange (coeffs).cast (x1_t {}), std::forward<Ts> (args)...);
    }
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  void reset_coeffs_on_idx (uint idx, Ts&&... args)
  {
    assert (idx < n_elems);
    get_target_coeffs<Part> (_coeffs[idx], std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part>
  void reset_states_on_idx (uint idx)
  {
    static_assert (mp11::mp_find<parts, Part>::value < n_parts, "");
    assert (idx < n_elems);
    Part::template reset_states<value_type> (_states[idx]);
  }
  //----------------------------------------------------------------------------
  // generic 0 setting, for safety use the variant above
  void reset_states_on_idx (uint idx)
  {
    assert (idx < n_elems);
    crange_memset (make_crange (_states[idx]), 0);
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    static_assert (mp11::mp_find<parts, Part>::value < n_parts, "");
    assert (idx < n_elems);
    return Part::template tick<value_type> (
      _coeffs[idx], _states[idx], std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, class... Ts>
  void reset_coeffs (Ts&&... args)
  {
    reset_coeffs_on_idx<Part> (0, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part>
  void reset_states()
  {
    reset_states_on_idx<Part> (0);
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
    static_assert (Idx < n_elems);
    reset_coeffs_on_idx<Part> (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <class Part, uint Idx>
  void reset_states()
  {
    static_assert (Idx < n_elems);
    reset_states_on_idx<Part> (Idx);
  }
  //----------------------------------------------------------------------------
  template <class Part, uint Idx, class... Ts>
  auto tick (Ts&&... args)
  {
    static_assert (Idx < n_elems);
    return tick_on_idx<Part> (Idx, std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  auto& get_coeffs()
  {
    static_assert (Idx < n_elems, "");
    return _coeffs[Idx];
  }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<coeff_type> get_coeffs (uint idx = 0)
  {
    assert (idx < n_elems);
    return (idx < n_elems) ? make_crange (_coeffs[idx]) : crange<value_type> {};
  }
  //----------------------------------------------------------------------------
  // for bulk coefficient smoothing
  crange<coeff_type> get_all_coeffs()
  {
    return {_coeffs.data(), _coeffs.size() * _coeffs[0].size()};
  }
  //----------------------------------------------------------------------------
  void reset_states() { crange_memset (make_crange (_states), 0); }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  using to_n_coeffs = k_int<T::n_coeffs>;
  template <class T>
  using to_n_states = k_int<T::n_states>;

  using n_coeffs_types = mp11::mp_transform<to_n_coeffs, parts>;
  using n_states_types = mp11::mp_transform<to_n_states, parts>;

  static constexpr uint max_coeffs
    = mp11::mp_max_element<n_coeffs_types, mp11::mp_less>::value;

  static constexpr uint max_states
    = mp11::mp_max_element<n_states_types, mp11::mp_less>::value;

  alignas (sizeof (value_type))
    std::array<std::array<coeff_type, max_coeffs>, n_elems> _coeffs;
  alignas (sizeof (value_type))
    std::array<std::array<value_type, max_states>, n_elems> _states;
};
} // namespace detail

//------------------------------------------------------------------------------
// N instances of one of many different DSP parts. Every instance can be only
// one of the available DSP parts simultaneously. This is controlled externally.
template <class Vect, uint Count, class... Parts>
using parts_to_class_one_of
  = detail::parts_to_class_one_of<Vect, Vect, Count, Parts...>;
//------------------------------------------------------------------------------
// As "parts_to_class_one_of" but the coefficient set is of width = 1,
// independently of the with of "Vect"
template <class Vect, uint Count, class... Parts>
using parts_to_class_one_of_single_coeff = detail::
  parts_to_class_one_of<Vect, vec_value_type_t<Vect>, Count, Parts...>;
//------------------------------------------------------------------------------
} // namespace artv
