#pragma once

#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
//------------------------------------------------------------------------------
struct no_interp {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 0;
  static constexpr uint n_coeffs_int = 0;
  static constexpr uint n_states     = 0;
  static constexpr uint n_points     = 1;
  static constexpr uint x_offset     = 0;
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_coeffs (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (std::array<V, n_points> y, V x)
  {
    return y[0];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, std::array<V, n_points> y, V x)
  {
    return tick (y, x);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct linear_interp {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 0;
  static constexpr uint n_coeffs_int = 0;
  static constexpr uint n_states     = 0;
  static constexpr uint n_points     = 2;
  static constexpr uint x_offset     = 0;
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_coeffs (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (std::array<V, n_points> y, V x)
  {
    return y[0] + x * (y[1] - y[0]);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, std::array<V, n_points> y, V x)
  {
    return tick (y, x);
  }
  //----------------------------------------------------------------------------
};

template <uint N>
struct lagrange_interp;

//------------------------------------------------------------------------------
template <>
struct lagrange_interp<2> {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 0;
  static constexpr uint n_coeffs_int = 0;
  static constexpr uint n_states     = 0;
  static constexpr uint n_points     = 3;
  static constexpr uint x_offset     = 0;
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_coeffs (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (std::array<V, n_points> y, V x)
  {
    using T = vec_value_type_t<V>;
    // terms with x

    V x0 = x;
    V x1 = x - vec_set<V> ((T) 1.);
    V x2 = x - vec_set<V> ((T) 2.);

    // A modern optimizer sees if factoring or parallelizing is better. Kept
    // naive.
    auto c1 = x1 * x2 * vec_set<V> ((T) 0.5);
    auto c2 = -(x0 * x2);
    auto c3 = x0 * x1 * vec_set<V> ((T) 0.5);
    return y[0] * c1 + y[1] * c2 + y[2] * c3;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, std::array<V, n_points> y, V x)
  {
    return tick (y, x);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <>
struct lagrange_interp<3> {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 0;
  static constexpr uint n_coeffs_int = 0;
  static constexpr uint n_states     = 0;
  static constexpr uint n_points     = 4;
  static constexpr uint x_offset     = 0;
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_coeffs (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (std::array<V, n_points> y, V x)
  {
    using T = vec_value_type_t<V>;

    V x0 = x;
    V x1 = x - vec_set<V> ((T) 1.);
    V x2 = x - vec_set<V> ((T) 2.);
    V x3 = x - vec_set<V> ((T) 3.);

    // A modern optimizer sees if factoring or parallelizing is better. Kept
    // naive.
    auto c1 = x1 * x2 * x3 * (-1.f / 6.f);
    auto c2 = x0 * x2 * x3 * 0.5f;
    auto c3 = x0 * x1 * x3 * -0.5f;
    auto c4 = x0 * x1 * x2 * (1. / 6.f);
    return y[0] * c1 + y[1] * c2 + y[2] * c3 + y[3] * c4;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, std::array<V, n_points> y, V x)
  {
    return tick (y, x);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct hermite_interp {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 0;
  static constexpr uint n_coeffs_int = 0;
  static constexpr uint n_states     = 0;
  static constexpr uint n_points     = 4;
  static constexpr uint x_offset     = 1;
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_coeffs (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (std::array<V, n_points> y, V x)
  {
    using T = vec_value_type_t<V>;

    V a = (((T) 3 * (y[1] - y[2])) - y[0] + y[3]) * (T) 0.5;
    V b = y[2] + y[2] + y[0] - ((T) 5 * y[1] + y[3]) * (T) 0.5;
    V c = (y[2] - y[0]) * (T) 0.5;
    return ((a * x + b) * x + c) * x + y[1];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, std::array<V, n_points> y, V x)
  {
    return tick (y, x);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct catmull_rom_interp {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs     = 0;
  static constexpr uint n_coeffs_int = 0;
  static constexpr uint n_states     = 0;
  static constexpr uint n_points     = 4;
  static constexpr uint x_offset     = 1; // interpolates between y[1] and y[2]
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_coeffs (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (std::array<V, n_points> y, V x)
  {
    using T = vec_value_type_t<V>;

    V x2 = x * x;
    V x3 = x2 * x;

    V r = x3 * (-y[0] + (T) 3 * y[1] - (T) 3 * y[2] + y[3]) * (T) 0.5;
    r += x2 * ((T) 2 * y[0] - (T) 5 * y[1] + (T) 4 * y[2] - y[3]) * (T) 0.5;
    r += x * (-y[0] + y[2]) * (T) 0.5;
    r += y[1];

    return r;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, std::array<V, n_points> y, V x)
  {
    return tick (y, x);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

} // namespace artv
