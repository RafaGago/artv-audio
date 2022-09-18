#pragma once

#include <cassert>
#include <cstring>
#include <type_traits>

#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
// bulk SIMD smoothing of many parameters of the same floating point type via a
// 1 pole lowpass.
//
// T = arithmetic type
// T_struct, struct only containing members of the T type (unverifiable)
//
// Notice that this was previously not using structs, but indexes. It wasn't
// very convenient, so this was deemed as the lesser evil.

template <class T, class T_struct, uint Simd_size = sse_bytes>
class value_smoother {
public:
  //----------------------------------------------------------------------------
  static_assert (std::is_floating_point_v<T>);
  static_assert (std::is_standard_layout_v<T_struct>);
  static_assert (!std::is_arithmetic_v<T_struct>);
  // That the struct contains only T's can't be validated, just trying...
  static_assert ((sizeof (T_struct) % sizeof (T)) == 0);
  //----------------------------------------------------------------------------
  using value_type                = T_struct;
  using arith_type                = T;
  static constexpr uint simd_size = Simd_size;
  static constexpr uint size      = sizeof (T_struct) / sizeof (T);
  using vec_type                  = vec<T, simd_size / sizeof (T)>;
  static constexpr uint vec_size  = vec_traits_t<vec_type>::size;
  //----------------------------------------------------------------------------
  void reset (T t_spl, T hz = 10.f)
  {
    memset (&_target, 0, sizeof _target);
    memset (&_current, 0, sizeof _current);
    using x1_t = vec<T, 1>;
    onepole_smoother::reset_coeffs (
      make_xspan (_smoother_coeff).cast (x1_t {}), vec_set<x1_t> (hz), t_spl);
  }
  //----------------------------------------------------------------------------
  void set_all_from_target() { _current = _target; }
  //----------------------------------------------------------------------------
  value_type const& get() const { return _current.value; }
  //----------------------------------------------------------------------------
  value_type& target() { return _target.value; }
  //----------------------------------------------------------------------------
  void tick (uint samples = 1)
  {
    for (uint i = 0; i < samples; ++i) {
      for (uint j = 0; j < _target.mem.size(); ++j) {
        _current.mem[j] = onepole_smoother::tick<vec_type> (
          make_xspan (_smoother_coeff),
          make_xspan (_current.mem[j]),
          _target.mem[j]);
      }
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  union values_union {
    value_type value;
    alignas (sse_bytes) simd_vec_array<T, size, sse_bytes> mem;
  };
  //----------------------------------------------------------------------------
  values_union _target;
  values_union _current;
  T            _smoother_coeff;
};
//------------------------------------------------------------------------------

} // namespace artv
