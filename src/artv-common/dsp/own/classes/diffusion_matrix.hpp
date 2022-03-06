#pragma once

#include <array>
#include <gcem.hpp>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <uint N>
struct householder_matrix {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);

    using T = vec_value_type_t<V>;
    std::array<V, N> y;

    V factor = vec_set<V> ((T) 0);
    for (uint i = 0; i < N; ++i) {
      factor += x[i];
    }
    factor *= ((T) (2. / (double) N));
    for (uint i = 0; i < N; ++i) {
      y[i] = (x[i] - factor);
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// https://signalsmith-audio.co.uk/writing/2021/lets-write-a-reverb/
// Gist there:
// https://gist.github.com/geraintluff/c55f9482dce0db6e6a1f4509129e9a2a
// I like it how much more clear it becomes VS the Wikipedia algorithm when made
// recursive.
template <uint N>
class hadamard_matrix {
public:
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);
    using T                      = vec_value_type_t<V>;
    static constexpr double norm = 1. / gcem::sqrt ((double) N);

    auto y = tick_raw (x);
    for (auto& v : y) {
      v *= (T) norm;
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <uint>
  friend class hadamard_matrix;

  template <class V>
  static std::array<V, N> tick_raw (crange<const V> x)
  {
    if constexpr (N > 1) {
      constexpr int half_n = N / 2;
      // all this superfluous array copying stuff should be removed by the
      // optimizer.
      auto y = array_join (
        hadamard_matrix<half_n>::template tick_raw<V> (x.get_head (half_n)),
        hadamard_matrix<half_n>::template tick_raw<V> (x.advanced (half_n)));

      for (int i = 0; i < half_n; ++i) {
        V a           = y[i];
        V b           = y[i + half_n];
        y[i]          = (a + b);
        y[i + half_n] = (a - b);
      }
      return y;
    }
    else {
      return std::array<V, 1> {x[0]};
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint N>
struct almost_hadamard_matrix;

// source Thierry Rochebois:
// https://github.com/axoloti/axoloti-contrib/blob/1.0.12/objects/tiar/FDN/AH6.axo
template <>
struct almost_hadamard_matrix<6> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 6;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);
    static constexpr double norm = 1. / 3. * gcem::sqrt (2.);

    using T = vec_value_type_t<V>;
    std::array<V, N> y;

    auto yA = x[0] + x[1] + x[2];
    auto yB = x[3] + x[4] + x[5];

    auto a = yA - x[0];
    auto b = yB - x[3];

    y[0] = a + b;
    y[3] = a - b;

    a    = yA - x[1];
    b    = yB - x[4];
    y[1] = a + b;
    y[4] = a - b;

    a    = yA - x[2];
    b    = yB - x[5];
    y[2] = a + b;
    y[5] = a - b;

    for (uint i = 0; i < N; ++i) {
      y[i] *= (T) norm;
    }

    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};

// source Thierry Rochebois:
// https://github.com/axoloti/axoloti-contrib/blob/1.0.12/objects/tiar/FDN/AH7.axo
template <>
struct almost_hadamard_matrix<7> {
  static constexpr uint N = 7;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);
    static constexpr double A = -0.26120387496374144251476820691706;
    static constexpr double B = 0.4459029062228060818860761551878;

    using T = vec_value_type_t<V>;
    std::array<V, N> y;

    y[0] = (x[1] + x[3] + x[5]) * (T) A + (x[0] + x[2] + x[4] + x[6]) * (T) B;
    y[1] = (x[0] + x[3] + x[4]) * (T) A + (x[1] + x[2] + x[5] + x[6]) * (T) B;
    y[2] = (x[2] + x[3] + x[6]) * (T) A + (x[0] + x[1] + x[4] + x[5]) * (T) B;
    y[3] = (x[0] + x[1] + x[2]) * (T) A + (x[3] + x[4] + x[5] + x[6]) * (T) B;
    y[4] = (x[1] + x[4] + x[6]) * (T) A + (x[0] + x[2] + x[3] + x[5]) * (T) B;
    y[5] = (x[0] + x[5] + x[6]) * (T) A + (x[1] + x[2] + x[3] + x[4]) * (T) B;
    y[6] = (x[2] + x[4] + x[5]) * (T) A + (x[0] + x[1] + x[3] + x[6]) * (T) B;

    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint N>
struct rotation_matrix;

template <>
struct rotation_matrix<4> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 4;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x, crange<std::array<V, 2>> w)
  {
    assert (x.size() >= N);
    assert (w.size() >= N / 4);

    //     1 2 3 4
    // ------------
    // 1 | + - - +
    // 2 | + - + -
    // 3 | + + - -
    // 4 | + + + +

    std::array<V, N> y;

    auto w1 = w[0][0];
    auto w2 = w[0][1];

    auto y1 = w1 * x[0] - w2 * x[1];
    auto y2 = w2 * x[0] + w1 * x[1];

    auto y3 = w1 * x[2] - w2 * x[3];
    auto y4 = w2 * x[2] + w1 * x[3];

    y[0] = w1 * y1 - w2 * y3;
    y[1] = w2 * y1 + w1 * y3;

    y[2] = w1 * y2 - w2 * y4;
    y[3] = w2 * y2 + w1 * y4;

    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};

template <>
struct rotation_matrix<8> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 8;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x, crange<std::array<V, 2>> w)
  {
    assert (x.size() >= N);
    assert (w.size() >= N / 4);

    //    1 2 3 4
    // ------------
    // 1 | + - - +
    // 2 | + - + -
    // 3 | + + - -
    // 4 | + + + +

    //     5 6 7 8
    // ------------
    // 5 | + - - +
    // 6 | + - + -
    // 7 | + + - -
    // 8 | + + + +

    //     1 2 3 4 5 6 7 8
    // ---------------------
    // 1 | + - - + - + + -
    // 2 | + - - + + - - +
    // 3 | + - + - - + - +
    // 4 | + - + - + - + -
    // 5 | + + - - - - + +
    // 6 | + + - - + + - -
    // 7 | + + + + - - - -
    // 8 | + + + + + + + +

    std::array<V, N> y;

    auto y_a = rot_matrix_4 (x.get_head (4), w.get_head (1));
    auto y_b = rot_matrix_4 (x.advanced (4), w.advanced (1));

    auto y1 = y_a[0];
    auto y2 = y_a[1];
    auto y3 = y_a[2];
    auto y4 = y_a[3];
    auto y5 = y_b[0];
    auto y6 = y_b[1];
    auto y7 = y_b[2];
    auto y8 = y_b[3];

    auto w1a = w[0][0];
    auto w2a = w[0][1];
    auto w1b = w[1][0];
    auto w2b = w[1][1];

    y[0] = w1a * y1 - w2a * y5;
    y[1] = w2a * y1 + w1a * y5;
    y[2] = w1a * y2 - w2a * y6;
    y[3] = w2a * y2 + w1a * y6;
    y[4] = w1b * y3 - w2b * y7;
    y[5] = w2b * y3 + w1b * y7;
    y[6] = w1b * y4 - w2b * y8;
    y[7] = w2b * y4 + w1b * y8;

    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};

template <>
struct rotation_matrix<16> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 16;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x, crange<std::array<V, 2>> w)
  {
    assert (x.size() >= N);
    assert (w.size() >= N / 4);

    //     1 2 3 4 5 6 7 8
    // ---------------------
    // 1 | + - - + - + + -
    // 2 | + - - + + - - +
    // 3 | + - + - - + - +
    // 4 | + - + - + - + -
    // 5 | + + - - - - + +
    // 6 | + + - - + + - -
    // 7 | + + + + - - - -
    // 8 | + + + + + + + +

    //      9 0 1 2 3 4 5 6
    // ---------------------
    // 9  | + - - + - + + -
    // 10 | + - - + + - - +
    // 11 | + - + - - + - +
    // 12 | + - + - + - + -
    // 13 | + + - - - - + +
    // 14 | + + - - + + - -
    // 15 | + + + + - - - -
    // 16 | + + + + + + + +

    //      1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
    //    ---------------------------------
    //  1 | + - - + - + + - - + + - + - - +
    //  2 | + - - + - + + - + - - + - + + -
    //  3 | + - - + + - - + - + + - - + + -
    //  4 | + - - + + - - + + - - + + - - +
    //  5 | + - + - - + - + - + - + + - + -
    //  6 | + - + - - + - + + - + - - + - +
    //  7 | + - + - + - + - - + - + - + - +
    //  8 | + - + - + - + - + - + - + - + -
    //  9 | + + - - - - + + - - + + + + - -
    // 10 | + + - - - - + + + + - - - - + +
    // 11 | + + - - + + - - - - + + - - + +
    // 12 | + + - - + + - - + + - - + + - -
    // 13 | + + + + - - - - - - - - + + + +
    // 14 | + + + + - - - - + + + + - - - -
    // 15 | + + + + + + + + - - - - - - - -
    // 16 | + + + + + + + + + + + + + + + +

    std::array<V, N> y;

    auto y_a = rot_matrix_8 (x.get_head (8), w.get_head (2));
    auto y_b = rot_matrix_8 (x.advanced (8), w.advanced (2));

    auto y1 = y_a[0];
    auto y2 = y_a[1];
    auto y3 = y_a[2];
    auto y4 = y_a[3];
    auto y5 = y_a[4];
    auto y6 = y_a[5];
    auto y7 = y_a[6];
    auto y8 = y_a[7];
    auto y9 = y_b[0];
    auto ya = y_b[1];
    auto yb = y_b[2];
    auto yc = y_b[3];
    auto yd = y_b[4];
    auto ye = y_b[5];
    auto yf = y_b[6];
    auto yg = y_b[7];

    auto w1a = w[0][0];
    auto w2a = w[0][1];
    auto w1b = w[1][0];
    auto w2b = w[1][1];
    auto w1c = w[2][0];
    auto w2c = w[2][1];
    auto w1d = w[3][0];
    auto w2d = w[3][1];

    y[0]  = w1a * y1 - w2a * y9;
    y[1]  = w2a * y1 + w1a * y9;
    y[2]  = w1a * y2 - w2a * ya;
    y[3]  = w2a * y2 + w1a * ya;
    y[4]  = w1b * y3 - w2b * yb;
    y[5]  = w2b * y3 + w1b * yb;
    y[6]  = w1b * y4 - w2b * yc;
    y[7]  = w2b * y4 + w1b * yc;
    y[8]  = w1c * y5 - w2c * yd;
    y[9]  = w2c * y5 + w1c * yd;
    y[10] = w1c * y6 - w2c * ye;
    y[11] = w2c * y6 + w1c * ye;
    y[12] = w1d * y7 - w2d * yf;
    y[13] = w2d * y7 + w1d * yf;
    y[14] = w1d * y8 - w2d * yg;
    y[15] = w2d * y8 + w1d * yg;

    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// TODO; there are more experimental matrices on "_reverb.jsfx-inc"
//------------------------------------------------------------------------------
// "rochebois_matrix" are conference matrices. Some experimentation was done
// with of them are found on "_reverb.jsfx-inc" but they are not ported. Albeit
// more expensive, the rotation matrices sounded better for me.
//------------------------------------------------------------------------------
template <uint N>
struct rochebois_matrix;

// source Thierry Rochebois:
// https://github.com/axoloti/axoloti-contrib/blob/1.0.12/objects/tiar/FDN/D4.axo
template <>
struct rochebois_matrix<4> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 4;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);
    static constexpr double norm = 1. / gcem::sqrt (3.);

    using T = vec_value_type_t<V>;
    std::array<V, N> y;

    y[0] = x[1] + x[2] + x[3];
    y[1] = x[0] + x[2] - x[3];
    y[2] = x[0] - x[1] + x[3];
    y[3] = x[0] + x[1] - x[2];

    for (uint i = 0; i < N; ++i) {
      y[i] *= (T) norm;
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};

// source Thierry Rochebois:
// https://github.com/axoloti/axoloti-contrib/blob/1.0.12/objects/tiar/FDN/D6.axo
template <>
struct rochebois_matrix<6> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 6;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);
    static constexpr double norm = 1. / gcem::sqrt (5.);

    using T = vec_value_type_t<V>;
    std::array<V, N> y;

    auto tA = x[4] + x[5];
    auto tB = x[2] + x[3];
    y[0]    = tA + tB + x[1];
    y[1]    = tA - tB + x[0];

    tA   = x[0] - x[1];
    tB   = x[5] - x[4];
    y[2] = tA + tB + x[3];
    y[3] = tA - tB + x[2];

    tA   = x[0] + x[1];
    tB   = x[3] - x[2];
    y[4] = tA + tB - x[5];
    y[5] = tA - tB - x[4];

    for (uint i = 0; i < N; ++i) {
      y[i] *= norm;
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};

// source Thierry Rochebois:
// https://github.com/axoloti/axoloti-contrib/blob/1.0.12/objects/tiar/FDN/D8.axo
template <>
struct rochebois_matrix<8> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 8;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);
    static constexpr double norm = 1. / gcem::sqrt (7.);

    using T = vec_value_type_t<V>;
    std::array<V, N> y;

    auto t = x[3] + x[4] + x[5] + x[6] + x[7];

    y[0] = t - (x[4] * (T) 2) + x[2] - x[1];
    y[1] = t - (x[5] * (T) 2) + x[0] - x[2];
    y[2] = t - (x[6] * (T) 2) + x[1] - x[0];

    auto t2 = -x[0] - x[1] - x[2] - x[3];
    y[3]    = (t - (x[7] * (T) 2) + t2);

    y[4] = t2 + (x[0] * (T) 2) + x[5] - x[6] + x[7];
    y[5] = t2 + (x[1] * (T) 2) - x[4] + x[6] + x[7];
    y[6] = t2 + (x[2] * (T) 2) + x[4] - x[5] + x[7];
    y[7] = t2 + (x[3] * (T) 2) - x[4] - x[5] - x[6];

    for (uint i = 0; i < N; ++i) {
      y[i] *= (T) norm;
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};

// source Thierry Rochebois:
// https://github.com/axoloti/axoloti-contrib/blob/1.0.12/objects/tiar/FDN/D10.axo
template <>
struct rochebois_matrix<10> {
  //----------------------------------------------------------------------------
  static constexpr uint N = 10;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (crange<const V> x)
  {
    assert (x.size() >= N);
    static constexpr double norm = 1. / 3.; // 1/sqrt(9)

    using T = vec_value_type_t<V>;
    std::array<V, N> y;

    auto t  = x[5] + x[6] + x[7] + x[8] + x[9];
    y[0]    = t + x[3] + x[4] - x[1] - x[2] - (x[5] * (T) 2);
    auto t2 = x[3] - x[4];
    auto t3 = t - x[0];
    y[1]    = t3 + t2 + x[2] - (x[6] * (T) 2);
    y[2]    = t3 - t2 + x[1] - (x[7] * (T) 2);
    t2      = x[1] - x[2];
    t3      = t + x[0];
    y[3]    = t3 + t2 - x[4] - (x[8] * (T) 2);
    y[4]    = t3 - t2 - x[3] - (x[9] * (T) 2);

    t    = x[0] + x[1] + x[2] + x[3] + x[4];
    t2   = x[7] - x[9];
    t3   = t - x[8];
    y[5] = t3 + t2 + x[6] - (x[0] * (T) 2);
    y[6] = t3 - t2 + x[5] - (x[1] * (T) 2);
    t2   = x[5] - x[9];
    t3   = t - x[6];

    y[9] = t + x[6] + x[8] - x[5] - x[7] - (x[4] * (T) 2);
    y[7] = t3 + t2 + x[8] - (x[2] * (T) 2);
    y[8] = t3 - t2 + x[7] - (x[3] * (T) 2);

    for (uint i = 0; i < N; ++i) {
      y[i] *= (T) norm;
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static std::array<T, N> tick (crange<const T> x)
  {
    using vectype = vec<T, 1>;
    return vec_array_unwrap<1, T> (
      tick<vectype> (x.template cast<const vectype>()));
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
