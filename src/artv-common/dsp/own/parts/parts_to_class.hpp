#pragma once

#include <array>
#include <cassert>
#include <utility>

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
  template <class... Ts>
  void reset_coeffs_on_idx (uint idx, Ts&&... args)
  {
    assert (idx < n_elems);
    part::reset_coeffs (
      make_crange (_coeffs[idx]).cast (builtin {}), std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
  void reset_states_on_idx (uint idx)
  {
    assert (idx < n_elems);
    part::reset_states (make_crange (_states[idx]).cast (builtin {}));
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto tick_on_idx (uint idx, Ts&&... args)
  {
    assert (idx < n_elems);
    return part::tick (
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
    T out = in;
    for (uint i = 0; i < n_elems; ++i) {
      out = tick_on_idx (i, out, std::forward<Ts> (args)...);
    }
    return out;
  }
//----------------------------------------------------------------------------
#if 0
  // TODO: forward fix unsmoothable coeffs if external smoothing is required.
  //----------------------------------------------------------------------------
  std::array<value_type, part::n_coeffs>* get_coeffs (uint idx = 0)
  {
    assert (idx < n_elems);
    return (idx < n_elems) ? &_coeffs[idx] : nullptr;
  }
  //----------------------------------------------------------------------------
#endif
private:
  std::array<std::array<value_type, part::n_coeffs>, n_elems> _coeffs;
  std::array<std::array<value_type, part::n_states>, n_elems> _states;
};
//------------------------------------------------------------------------------

} // namespace artv
