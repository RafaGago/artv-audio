#pragma once

#include "artv-common/dsp/own/classes/windowed_sync.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {
//------------------------------------------------------------------------------
struct zero_order_hold {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs          = 0;
  static constexpr uint n_states          = 0;
  static constexpr uint n_points          = 1;
  static constexpr uint x_offset          = 0;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = true;
  static constexpr bool states_are_vec    = true;
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
  static constexpr uint n_coeffs          = 0;
  static constexpr uint n_states          = 0;
  static constexpr uint n_points          = 2;
  static constexpr uint x_offset          = 0;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = true;
  static constexpr bool states_are_vec    = true;
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
  static constexpr uint n_coeffs          = 0;
  static constexpr uint n_states          = 0;
  static constexpr uint n_points          = 3;
  static constexpr uint x_offset          = 0;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = true;
  static constexpr bool states_are_vec    = true;
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
  static constexpr uint n_coeffs          = 0;
  static constexpr uint n_states          = 0;
  static constexpr uint n_points          = 4;
  static constexpr uint x_offset          = 0;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = true;
  static constexpr bool states_are_vec    = true;
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
  static constexpr uint n_coeffs          = 0;
  static constexpr uint n_states          = 0;
  static constexpr uint n_points          = 4;
  static constexpr uint x_offset          = 1;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = true;
  static constexpr bool states_are_vec    = true;
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
  static constexpr uint n_coeffs = 0;
  static constexpr uint n_states = 0;
  static constexpr uint n_points = 4;
  static constexpr uint x_offset = 1; // interpolates between y[1] and y[2]
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = true;
  static constexpr bool states_are_vec    = true;
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
template <uint N, uint N_tables>
struct sinc_interp {
  //----------------------------------------------------------------------------
  static constexpr uint   n_coeffs          = N * (N_tables + 1);
  static constexpr uint   n_states          = 0;
  static constexpr uint   n_points          = N;
  static constexpr uint   x_offset          = N / 2;
  static constexpr double n_tables          = N_tables;
  static constexpr double mu                = 1. / n_tables;
  static constexpr bool   coeffs_are_vec    = false;
  static constexpr bool   coeffs_are_global = true;
  static constexpr bool   states_are_vec    = true;
  //----------------------------------------------------------------------------
  template <class T> // fc 0 to 0.5
  static void reset_coeffs (crange<T> co, float fc, float kaiser_att_db)
  {
    static_assert (!is_vec_v<T>);
    assert (co.size() >= n_coeffs);

    constexpr auto k_mu = N % 2 ? 0. : 0.5;

    auto mem = co;
    for (uint tbl = 0; tbl < (N_tables + 1); ++tbl) {
      auto tblm = mem.cut_head (n_points);
      kaiser_lp_kernel (
        tblm, (T) fc, (T) kaiser_att_db, (T) (-mu * tbl + k_mu), false);

      constexpr bool tbl_debug = false;
      if constexpr (tbl_debug) {
        printf ("tbl %u: ", tbl);
        for (float v : tblm) {
          printf ("%f ", v);
        }
        puts ("");
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co,
    crange<V>,
    std::array<V, n_points> y,
    V                       x)
  {
    using T = vec_value_type_t<V>;
    for (uint i = 1; i < vec_traits_t<V>::size; ++i) {
      // this interpolation method doesn't support multiple lookup points.
      assert (x[0] == x[i]);
    }
    // there are other implentations that use more memory at the expense of less
    // calculations. This is the "naive" way, going for more locality and half
    // the table size at the expense of more muls and additions (almost free on
    // a modern CPU). For float and N up to 8 it only touches one cache line.
    T    frac     = x[0];
    uint tbl      = frac * (n_tables);
    T    co_lerp2 = (frac - ((T) tbl * mu)) * n_tables;
    T    co_lerp1 = 1.f - co_lerp2;

    auto dotprod = vec_set<V> (0);
    for (uint i = 0; i < n_points; ++i) {
      // table interpolation
      T sinc_lerp = co[tbl * n_points + i] * co_lerp1;
      sinc_lerp += co[(tbl + 1) * n_points + i] * co_lerp2;
      dotprod += y[i] * sinc_lerp;
    }
    return dotprod;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
