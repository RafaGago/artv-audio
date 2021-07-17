#pragma once

/* Fiters from: https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf*/

#include <cmath>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace andy {

// TODO: add the version with better precission which uses "sin" for coeff
// calculation ?

struct svf {
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, m0, m1, m2, n_coeffs };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void lowpass (
    crange<T> c,
    double    freq,
    double    q,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 0.;
    *co.m1   = 0.;
    *co.m2   = 1.;
  }
  //----------------------------------------------------------------------------
  template <size_t simd_bytes, class T>
  static void lowpass_multi_aligned (
    crange<T>       c,
    crange<const T> freq, // no alignment required
    crange<const T> q, // no alignment required
    T               sr)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= (n_coeffs * n_builtins));
    assert (freq.size() >= n_builtins);
    assert (q.size() >= n_builtins);

    simdreg f {freq.data(), xsimd::unaligned_mode {}};
    simdreg g = xsimd::tan (f * ((T) M_PI) / sr);

    simdreg qsimd {q.data(), xsimd::unaligned_mode {}};
    simdreg k {((T) 1.0) / qsimd};

    simdreg a1_v {(T) 1.0};
    a1_v /= g * (g + k) + ((T) 1.);
    a1_v.store_aligned (&c[a1 * n_builtins]);

    simdreg a2_v = g * a1_v;
    a2_v.store_aligned (&c[a2 * n_builtins]);

    simdreg a3_v = g * a2_v;
    a3_v.store_aligned (&c[a3 * n_builtins]);

    simdreg m0_v {(T) 0.};
    m0_v.store_aligned (&c[m0 * n_builtins]);

    simdreg m1_v {(T) 0.};
    m1_v.store_aligned (&c[m1 * n_builtins]);

    simdreg m2_v {(T) 1.};
    m2_v.store_aligned (&c[m2 * n_builtins]);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void bandpass (
    crange<T> c,
    double    freq,
    double    q,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 0.;
    *co.m1   = 1.;
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void highpass (
    crange<T> c,
    double    freq,
    double    q,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -k;
    *co.m2   = -1.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void notch (
    crange<T> c,
    double    freq,
    double    q,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -k;
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void peak (
    crange<T> c,
    double    freq,
    double    q,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -k;
    *co.m2   = -2.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void allpass (
    crange<T> c,
    double    freq,
    double    q,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -2. * k;
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <size_t simd_bytes, class T>
  static void allpass_multi_aligned (
    crange<T>       c,
    crange<const T> freq, // no alignment required
    crange<const T> q, // no alignment required
    T               sr)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= (n_coeffs * n_builtins));
    assert (freq.size() >= n_builtins);
    assert (q.size() >= n_builtins);

    simdreg f {freq.data(), xsimd::unaligned_mode {}};
    simdreg g = xsimd::tan (f * ((T) M_PI) / sr);

    simdreg qsimd {q.data(), xsimd::unaligned_mode {}};
    simdreg k {((T) 1.0) / qsimd};

    simdreg a1_v {(T) 1.0};
    a1_v /= g * (g + k) + ((T) 1.);
    a1_v.store_aligned (&c[a1 * n_builtins]);

    simdreg a2_v = g * a1_v;
    a2_v.store_aligned (&c[a2 * n_builtins]);

    simdreg a3_v = g * a2_v;
    a3_v.store_aligned (&c[a3 * n_builtins]);

    simdreg m0_v {(T) 1.};
    m0_v.store_aligned (&c[m0 * n_builtins]);

    simdreg m1_v {(T) -2.};
    m1_v *= k;
    m1_v.store_aligned (&c[m1 * n_builtins]);

    simdreg m2_v {(T) 0.};
    m2_v.store_aligned (&c[m2 * n_builtins]);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void bell (
    crange<T> c,
    double    freq,
    double    q,
    double    belldB,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double A = pow (10.f, belldB / 40.f);
    double g = tan (M_PI * freq / sr);
    double k = 1.0 / (q * A);
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = k * (A * A - 1.);
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void low_shelf (
    crange<T> c,
    double    freq,
    double    q,
    double    belldB,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double A = pow (10.f, belldB / 40.f);
    double g = tan (M_PI * freq / sr) / sqrt (A);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = k * (A - 1.);
    *co.m2   = (A * A) - 1.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void high_shelf (
    crange<T> c,
    double    freq,
    double    q,
    double    belldB,
    double    sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    c.shrink_head (coeff_offset);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double A = pow (10.f, belldB / 40.f);
    double g = tan (M_PI * freq / sr) * sqrt (A);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = A * A;
    *co.m1   = k * (1. - A) * A;
    *co.m2   = 1. - (A * A);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {
    static_assert (std::is_floating_point<T>::value, "");
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T> s)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               v0)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (s.size() >= n_states);
    assert (c.size() >= n_coeffs);

    return calc (
      s[ic1eq], s[ic2eq], c[m0], c[m1], c[m2], c[a1], c[a2], c[a3], v0);
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, multiple outs. (E.g. stereo filter using double)
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick (
    crange<const T>                                      c, // coeffs
    std::array<crange<T>, simd_reg<T, simd_bytes>::size> s, // state
    std::array<T, simd_reg<T, simd_bytes>::size>         v0s)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= n_coeffs);
    for (uint i = 0; i < n_builtins; ++i) {
      assert (s[i].size() >= n_states);
      assert (s[i].size() >= n_states);
    }

    simdreg ic1eq_v, ic2eq_v, v0;

    for (uint i = 0; i < n_builtins; ++i) {
      ic1eq_v[i] = s[i][ic1eq];
      ic2eq_v[i] = s[i][ic2eq];
      v0[i]      = v0s[i];
    }

    simdreg out = calc<simdreg> (
      ic1eq_v,
      ic2eq_v,
      simdreg {c[m0]},
      simdreg {c[m1]},
      simdreg {c[m2]},
      simdreg {c[a1]},
      simdreg {c[a2]},
      simdreg {c[a3]},
      v0);

    for (uint i = 0; i < n_builtins; ++i) {
      s[i][ic1eq] = ic1eq_v[i];
      s[i][ic2eq] = ic2eq_v[i];
    }
    return out;
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double)
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_aligned (
    crange<const T> c, // coeffs
    crange<T>       s, // coefficients interleaved, ready to SIMD load
    crange<const T> v0s) // N inputs ready to SIMD load
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_builtins * n_states);
    assert (v0s.size() >= n_builtins);

    simdreg ic1eq_v, ic2eq_v, v0;

    ic1eq_v.load_aligned (&s[ic1eq * n_builtins]);
    ic2eq_v.load_aligned (&s[ic2eq * n_builtins]);
    v0.load_aligned (v0s.data());

    simdreg out = calc<simdreg> (
      ic1eq_v,
      ic2eq_v,
      simdreg {c[m0]},
      simdreg {c[m1]},
      simdreg {c[m2]},
      simdreg {c[a1]},
      simdreg {c[a2]},
      simdreg {c[a3]},
      v0);

    ic1eq_v.store_aligned (&s[ic1eq * n_builtins]);
    ic2eq_v.store_aligned (&s[ic2eq * n_builtins]);

    return out;
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T> c, // coeffs interleaved, ready to SIMD load
    crange<T>       s, // coefficients interleaved, ready to SIMD load
    crange<const T> v0s) // N inputs ready to SIMD load
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg                    = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = simdreg::size;

    assert (c.size() >= n_builtins * n_coeffs);
    assert (s.size() >= n_builtins * n_states);
    assert (v0s.size() >= n_builtins);

    simdreg m0_v, m1_v, m2_v, a1_v, a2_v, a3_v, ic1eq_v, ic2eq_v, v0;

    ic1eq_v.load_aligned (&s[ic1eq * n_builtins]);
    ic2eq_v.load_aligned (&s[ic2eq * n_builtins]);
    v0.load_aligned (v0s.data());

    m0_v.load_aligned (&c[m0 * n_builtins]);
    m1_v.load_aligned (&c[m1 * n_builtins]);
    m2_v.load_aligned (&c[m2 * n_builtins]);
    a1_v.load_aligned (&c[a1 * n_builtins]);
    a2_v.load_aligned (&c[a2 * n_builtins]);
    a3_v.load_aligned (&c[a3 * n_builtins]);

    simdreg out
      = svf::calc (ic1eq_v, ic2eq_v, m0_v, m1_v, m2_v, a1_v, a2_v, a3_v, v0);

    ic1eq_v.store_aligned (&s[ic1eq * n_builtins]);
    ic2eq_v.store_aligned (&s[ic2eq * n_builtins]);

    return out;
  }
  //----------------------------------------------------------------------------
private:
  template <class T>
  static T calc (
    T& ic1eq_v,
    T& ic2eq_v,
    T  m0_v,
    T  m1_v,
    T  m2_v,
    T  a1_v,
    T  a2_v,
    T  a3_v,
    T  v0)
  {
    T v3    = v0 - ic2eq_v;
    T v1    = a1_v * ic1eq_v + a2_v * v3;
    T v2    = ic2eq_v + a2_v * ic1eq_v + a3_v * v3;
    T out   = m0_v * v0 + m1_v * v1 + m2_v * v2;
    ic1eq_v = (v1 * 2.) - ic1eq_v;
    ic2eq_v = (v2 * 2.) - ic2eq_v;
    return out;
  }
  //----------------------------------------------------------------------------
  template <class T>
  struct coeff_ptrs {
    T* a1;
    T* a2;
    T* a3;
    T* m0;
    T* m1;
    T* m2;
  };
  //----------------------------------------------------------------------------
  template <class T>
  static coeff_ptrs<T> unpack_interleaved_coeffs (crange<T> c, uint separation)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (c.size() >= ((n_coeffs - 1) * separation) + 1);

    coeff_ptrs<T> ret;
    ret.a1 = &c[a1 * separation];
    ret.a2 = &c[a2 * separation];
    ret.a3 = &c[a3 * separation];
    ret.m0 = &c[m0 * separation];
    ret.m1 = &c[m1 * separation];
    ret.m2 = &c[m2 * separation];
    return ret;
  }
  //----------------------------------------------------------------------------
};
}} // namespace artv::andy
