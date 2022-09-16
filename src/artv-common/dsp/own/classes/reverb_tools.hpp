#pragma once

#include <cmath>
#include <gcem.hpp>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/primes.hpp"
#include "artv-common/misc/primes_table.hpp"
#include "artv-common/misc/random_table.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#include "artv-common/dsp/own/classes/delay_line.hpp"

namespace artv {

// Allpass with an externally managed multichannel delay line (see
// delay_line.hpp), so interpolation can be reused.
class allpass_fn {
public:
  //----------------------------------------------------------------------------
  template <class T>
  struct result {
    T out;
    T to_push;
  };
  //----------------------------------------------------------------------------
  template <class T>
  static result<T> tick (T in, T yn, T g_fb, T g_ffw)
  {
    T u = in + yn * g_fb;
    return {yn - u * g_ffw, u};
  }
  //----------------------------------------------------------------------------
  // index 0 is the result, 1 the value to push to a delay line, yn is the
  // delayed sample fetched from the queue
  template <class T>
  static result<T> tick (T in, T yn, T gain)
  {
    return tick (in, yn, gain, gain);
  }

  //----------------------------------------------------------------------------
  // tick all delay line channels in parallel
  template <class T, class Time_type, class Delay_line>
  static void tick (
    crange<T>               out,
    crange<const T>         in,
    crange<const Time_type> delay,
    crange<const T>         gain,
    Delay_line&             dl)
  {
    assert (out.size() >= dl.n_channels());
    assert (in.size() >= dl.n_channels());
    assert (delay.size() >= dl.n_channels());
    assert (gain.size() >= dl.n_channels());

    assert (dl.n_channels() <= 512); // limit VLA range
    T to_push[dl.n_channels()];

    for (uint i = 0; i < dl.n_channels(); ++i) {
      auto r     = tick (in[i], dl.get (delay[i], i), gain[i]);
      out[i]     = r.out;
      to_push[i] = r.to_push;
    }
    dl.push (make_crange (&to_push[0], dl.n_channels()));
  }
  //----------------------------------------------------------------------------
  // tick all delay line channels in serial
  template <class T, class Time_type, class Delay_line>
  static T tick (
    T                       in,
    crange<const Time_type> delay,
    crange<const T>         gain,
    Delay_line&             dl)
  {
    assert (delay.size() >= dl.n_channels());
    assert (gain.size() >= dl.n_channels());

    assert (dl.n_channels() <= 512); // limit VLA range
    T to_push[dl.n_channels()];

    T out = in;
    for (uint i = 0; i < dl.n_channels(); ++i) {
      auto r     = tick (out, dl.get (delay[i], i), gain[i]);
      out        = r.out;
      to_push[i] = r.to_push;
    }
    dl.push (make_crange (&to_push[0], dl.n_channels()));
    return out;
  }
  //----------------------------------------------------------------------------
  // special-casing the case of 1 channel to remove the loop
  template <class T, class Time_type, class Delay_line>
  static T tick (T in, Time_type delay, T gain, Delay_line& dl)
  {
    assert (dl.n_channels() == 1);
    auto r = tick (in, dl.get (delay, 0), gain);
    dl.push (make_crange (&r.to_push, 1));
    return r.out;
  }
  //----------------------------------------------------------------------------
};

namespace detail {

template <class V, class Circ_bufer>
class allpass {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem) { _mem.reset (mem); }
  //----------------------------------------------------------------------------
  V tick (V in, uint delay, V gain)
  {
    auto r = allpass_fn::tick<V> (in, _mem[delay], gain);
    _mem.push (r.to_push);
    return r.out;
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
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
class allpass_with_params {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem) { _ap.reset (mem); }
  //----------------------------------------------------------------------------
  void set (uint delay_samples, V gain)
  {
    set_time (delay_samples);
    set_gain (gain);
  }
  //----------------------------------------------------------------------------
  void set_time (uint samples) { _delay = samples; }
  //----------------------------------------------------------------------------
  void set_gain (V gain) { _gain = gain; }
  //----------------------------------------------------------------------------
  V tick (V in) { return _ap.tick (in, _delay, _gain); }
  //----------------------------------------------------------------------------
private:
  allpass<V, true> _ap;
  uint             _delay {};
  V                _gain {};
};
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static V delay_get_feedback_gain_for_time (
  vec_value_type_t<V> time_sec, // desired
  vec_value_type_t<V> gain, // db_to_gain (-60). for RT60
  vec_value_type_t<V> srate,
  V                   delay_spls) // delay line length
{
  using T = vec_value_type_t<V>;

  // n_cycles = (time_sec * srate) / delay_spls;
  // x^n_cycles = att_lin, x = pow(att_lin, 1/n_cycles)
  return vec_pow (gain, delay_spls / (time_sec * srate));
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static V delay_get_feedback_gain_for_rt60_time (
  vec_value_type_t<V> time_sec, // desired
  vec_value_type_t<V> srate,
  V                   delay_spls) // delay line length
{
  using T                          = vec_value_type_t<V>;
  static constexpr auto rt_60_gain = (T) constexpr_db_to_gain (-60.);
  static constexpr auto ln_rt_60   = (T) gcem::log (rt_60_gain);

  return vec_exp (ln_rt_60 * (delay_spls / (time_sec * srate)));
}

// a namespace class...
struct delay_length {
  //----------------------------------------------------------------------------
  static uint meters_to_samples (float m, float srate, bool ceil_round)
  {
    float v = (srate * m * (1.f / 343.f)); // 343 : sounspeed m/s
    return ceil_round ? (uint) std::ceil (v) : (uint) v;
  }
  //----------------------------------------------------------------------------
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
