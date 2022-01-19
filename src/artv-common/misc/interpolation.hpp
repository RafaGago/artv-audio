#pragma once

namespace artv {

//------------------------------------------------------------------------------
struct no_interp {
  static constexpr uint n_points = 1;

  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V get (std::array<V, n_points> y, V x)
  {
    return y[0];
  }
};
//------------------------------------------------------------------------------
struct linear_interp {
  static constexpr uint n_points = 2;

  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V get (std::array<V, n_points> y, V x)
  {
    using T = vec_value_type_t<V>;
    return y[0] * (vec_set<V> ((T) 1) - x) + y[1] * x;
  }
};
//------------------------------------------------------------------------------
// _esi: equally spaced on integer boundaries
struct lagrange_interp_2_esi {
  static constexpr uint n_points = 3;

  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V get (std::array<V, n_points> y, V x)
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
};
//------------------------------------------------------------------------------
// _esi: equally spaced on integer boundaries
struct lagrange_interp_3_esi {
  static constexpr uint n_points = 4;

  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V get (std::array<V, n_points> y, V x)
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
};
//------------------------------------------------------------------------------
// _esi: equally spaced on integer boundaries
struct cubic_interp_esi {
  static constexpr uint n_points = 4;

  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V get (std::array<V, n_points> y, V x)
  {
    using T = vec_value_type_t<V>;

    static constexpr T co[] = {
      (T) -0.5,
      (T) 1.0,
      (T) -0.5,
      (T) 0.0,
      (T) 1.5,
      (T) -2.5,
      (T) 0.0,
      (T) 1.0,
      (T) -1.5,
      (T) 2.0,
      (T) 0.5,
      (T) 0.0,
      (T) 0.5,
      (T) -0.5,
      (T) 0.0,
      (T) 0.0};

    V x3 = vec_set<V> ((T) 1.);
    V x2 = x;
    V x1 = x2 * x2;
    V x0 = x2 * x1;

    // A modern optimizer sees if factoring or parallelizing is better. Kept
    // naive.

    V c0 = co[0] * x0 + co[1] * x1 + co[2] * x2 + co[3] * x3;
    V c1 = co[4] * x0 + co[5] * x1 + co[6] * x2 + co[7] * x3;
    V c2 = co[8] * x0 + co[9] * x1 + co[10] * x2 + co[11] * x3;
    V c3 = co[12] * x0 + co[13] * x1 + co[14] * x2 + co[15] * x3;

    return c0 * y[0] + c1 * y[1] + c2 * y[2] + c3 * y[3];
  }
};
//------------------------------------------------------------------------------

} // namespace artv
