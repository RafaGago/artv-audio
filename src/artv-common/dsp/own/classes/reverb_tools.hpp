#pragma once

#include <cmath>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/primes.hpp"
#include "artv-common/misc/primes_table.hpp"
#include "artv-common/misc/random_table.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#include "artv-common/dsp/own/classes/delay_line.hpp"

namespace artv {

namespace detail {
template <
  class V,
  class Circ_bufer,
  enable_if_vec_of_float_point_t<V>* = nullptr>
class allpass {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem) { _mem.reset (mem); }
  //----------------------------------------------------------------------------
  V tick (V in, uint delay, V gain)
  {
    V yn = _mem[delay];
    V y  = in + yn * gain;
    _mem.push (y);
    return yn - y * gain;
  }
  //----------------------------------------------------------------------------
  V tick (V in, uint delay, V gain_forward, V gain_backward)
  {
    V yn = _mem[delay];
    V y  = in + yn * gain_backward;
    _mem.push (y);
    return yn - y * gain_forward;
  }
  //----------------------------------------------------------------------------
private:
  Circ_bufer _mem;
};

} // namespace detail
//------------------------------------------------------------------------------
template <class V, bool Use_pow2_opt = true>
using allpass = detail::allpass<
  V,
  std::
    conditional_t<Use_pow2_opt, pow2_circular_buffer<V>, circular_buffer<V>>>;
//------------------------------------------------------------------------------
// Basic reverb building block. Strictly not a delay line. Not sure if it
// belongs here. It probably needs to be moved.

//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static constexpr V delay_get_feedback_gain_for_time (
  vec_value_type_t<V> time_sec, // desired
  vec_value_type_t<V> att, // -60. for RT60
  vec_value_type_t<V> srate,
  V                   delays_spls)
{
  using T = vec_value_type_t<V>;

  auto rate = srate / delays_spls;
  return vec_exp ((T) M_LN10 * ((T) 1 / (T) 20) * att / (rate * time_sec));
}
// a namespace class...
struct delay_length {
  //----------------------------------------------------------------------------
  static uint meters_to_samples (float m, float srate, bool ceil_round)
  {
    float v = (srate * m * (1.f / 343.f)); // 343 : sounspeed m/s
    return ceil_round ? (uint) std::ceil (v) : (uint) v;
  }
  //------------------------------------------------------------------------------
  // Gets delay leghts as pure prime numbers
  // work mem should ideally contain "primes_table_size_guess (spls_min,
  // spls_max)" samples
  //
  // The function can fail if there are not "n_delays" prime numbers between
  // "spls_min" and "spls_max" or if "work_mem" can not allocate enough space
  // for the intermediate prime table.
  template <class T>
  static bool get_prime (
    crange<T> dst,
    T         spls_min,
    T         spls_max,
    crange<T> work_mem)
  {
    assert (dst.size());
    auto tbl = make_primes_table (work_mem, spls_min, spls_max);
    if (unlikely (!tbl)) {
      assert (false);
      return false;
    }
    if (tbl.size() < dst.size()) {
      return false;
    }
    auto dist = tbl.size() / (dst.size() - 1);
    auto ret  = work_mem.get_head (dst.size());
    uint i    = 0;
    for (; i < dst.size() - 1; ++i) {
      dst[i] = (T) tbl[i * dist];
    }
    dst[i] = (T) tbl.last();
    return true;
  }
  //----------------------------------------------------------------------------
  // https://ccrma.stanford.edu/~jos/pasp/Prime_Power_Delay_Line.html
  // This has prime_idx to get an offset to the primes list and a rounding
  // factor to allow interpolated delay lines. The conversion from meters to
  // samples is the free function "delay_length_meters_to_samples".
  //----------------------------------------------------------------------------
  template <class T>
  static void get_prime_power (
    crange<T> dst,
    T         spls_min,
    T         spls_max,
    uint      prime_idx,
    float     rounding_factor)
  {
    assert (dst.size());
    float ln_ratio = log ((float) spls_min / (float) spls_max);

    for (uint i = 0; i < dst.size(); ++i) {
      float f = exp (((float) i / (float) (dst.size() - 1)) * ln_ratio);
      float v = spls_min * f;
      v = (log (v) / log (primes_table[prime_idx + i]) * rounding_factor) + 0.5;
      v = floor (v);
      v /= rounding_factor;
      dst[i] = (T) pow (primes_table[prime_idx + i], v);
    }
  }
  //------------------------------------------------------------------------------
  // A "get_prime" with "get_prime_power" as a fallback
  template <class T>
  static void get (
    crange<T> dst,
    T         spls_min,
    T         spls_max,
    uint      prime_idx,
    float     rounding_factor,
    crange<T> work_mem)
  {
    assert (dst.size());
    if (get_prime (dst, spls_min, spls_max, work_mem)) {
      return;
    }
    get_prime_power (dst, spls_min, spls_max, prime_idx, rounding_factor);
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
