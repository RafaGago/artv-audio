#pragma once

#include <array>
#include <gcem.hpp>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/misc/interpolation.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <uint N>
struct householder_matrix {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
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
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <uint>
  friend class hadamard_matrix;

  template <class V>
  static std::array<V, N> tick_raw (std::array<V, N> const& x)
  {
    if constexpr (N > 1) {
      constexpr int half_n = N / 2;
      // all this superfluous array copying stuff should be removed by the
      // optimizer.
      auto y = array_join (
        hadamard_matrix<half_n>::tick_raw (array_slice<0, half_n> (x)),
        hadamard_matrix<half_n>::tick_raw (array_slice<half_n, half_n> (x)));

      for (int i = 0; i < half_n; ++i) {
        V a           = y[i];
        V b           = y[i + half_n];
        y[i]          = (a + b);
        y[i + half_n] = (a - b);
      }
      return y;
    }
    else {
      return x;
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
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
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
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
  }
  //----------------------------------------------------------------------------
};
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
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
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
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
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
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
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
  static std::array<V, N> tick (std::array<V, N> x)
  {
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
  static std::array<T, N> tick (std::array<T, N> x)
  {
    return vec_array_unwrap<1, T> (tick (vec_array_wrap<1> (x)));
  }
  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
} // namespace artv
