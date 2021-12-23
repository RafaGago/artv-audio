#pragma once

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
// All the classes under parts are coded stateless for $REASONS (see README).
// This is a wrapper to make most of them classes.
//
// As of now it can operate host many instance of the same class. Instances
// are selected by passing a std::integral_constant/k_uint, this was to avoid
// collisions with the variadic parameter pack. Might be revised depending on
// how it turns out to be.
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
      coeffs.cast (builtin {}), std::forward<Ts> (args)...);
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
    part::template reset_states<value_type> (
      make_crange (_states[idx]).cast (builtin {}));
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    assert (idx < n_elems);
    return part::template tick<value_type> (
      make_crange (_coeffs[idx]).cast (builtin {}),
      make_crange (_states[idx]).cast (builtin {}),
      std::forward<Ts> (args)...);
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
  std::array<std::array<value_type, part::n_coeffs>, n_elems> _coeffs;
  std::array<std::array<value_type, part::n_states>, n_elems> _states;
};
//------------------------------------------------------------------------------
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
      coeffs.cast (builtin {}), std::forward<Ts> (args)...);
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
    part::template reset_states<value_type> (
      make_crange (_states).cast (builtin {}));
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
      make_crange (_coeffs).cast (builtin {}),
      make_crange (_states).cast (builtin {}),
      std::forward<Ts> (args)...);
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

  std::tuple<coeffs_array<Parts>...> _coeffs;
  std::tuple<states_array<Parts>...> _states;
};
//------------------------------------------------------------------------------
} // namespace artv
