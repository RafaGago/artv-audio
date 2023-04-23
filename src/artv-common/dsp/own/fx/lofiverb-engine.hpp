#pragma once

// internals of lofiverb.hpp

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "artv-common/dsp/own/classes/noise.hpp"
#include "artv-common/dsp/own/fx/lofiverb-algorithms.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/float.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_approx_math.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/vec_util.hpp"
#include "artv-common/misc/xspan.hpp"
#include "boost/mp11/algorithm.hpp"

namespace artv { namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
// truncating fixed point type for computation
using fixpt_tt = fixpt_d<1, 0, 15, 0>;
// rounding fixed point type for computation
using fixpt_tr = fixpt_d<1, 0, 15, fixpt_rounding>;

using fixpt_t = fixpt_tt;

using fixpt_acum_t = fixpt_s<1, 0, 31>;

// fixed point type for storage
using fixpt_sto      = fixpt_s<1, 0, 15, 0>;
using fixpt_spls     = fixpt_m<0, 14, 0, 0>;
using fixpt_spls_mod = fixpt_m<0, 9, 0, 0>;

using float16 = f16pack<5, -1, f16pack_dftz | f16pack_clamp>;
//------------------------------------------------------------------------------
enum class interpolation : u8 { thiran, linear };

struct free_storage_data {
  uint count;
};

struct allpass_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g;
  interpolation  interp;
};

struct comb_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g;
  interpolation  interp;
};

// 1 tap delay
struct delay_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  interpolation  interp;
};

struct variable_delay_data {
  fixpt_spls min_spls;
  fixpt_spls max_spls;
  bool       mod_is_abs;
};

struct block_delay_data {
  fixpt_spls spls;
  fixpt_spls extra_spls; // In samples. To be able to do negative offsets
};

struct multitap_delay_data {
  std::array<fixpt_spls, 12> spls;
  uint                       count;
};

// A multitap delay summed in parallel
struct parallel_delay_data : public multitap_delay_data {
  std::array<fixpt_t, 12> g;
  uint                    n_outs;
};

struct filter_data {
  fixpt_t g;
  bool    is_lowpass;
};

struct crossover_data {
  static constexpr uint                max_n_bands = 2;
  std::array<fixpt_t, max_n_bands>     g; // fc
  std::array<fixpt_t, max_n_bands + 1> band_gain;
  uint                                 n_bands;
};

struct quantizer_data {
  fixpt_t g;
};

using stage_data = std::variant<
  allpass_data,
  comb_data,
  delay_data,
  block_delay_data,
  multitap_delay_data,
  parallel_delay_data,
  filter_data,
  crossover_data,
  quantizer_data,
  variable_delay_data,
  free_storage_data>;
//------------------------------------------------------------------------------
static constexpr stage_data make_ap (
  u16           spls,
  float         g      = 0.f,
  u16           mod    = 0,
  interpolation interp = interpolation::thiran)
{
  return allpass_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g),
    interp};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_comb (
  u16           spls,
  float         g      = 0.f,
  u16           mod    = 0,
  interpolation interp = interpolation::thiran)
{
  return comb_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g),
    interp};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_delay (
  u16           spls,
  u16           mod    = 0,
  interpolation interp = interpolation::thiran)
{
  return delay_data {
    fixpt_spls::from_int (spls), fixpt_spls_mod::from_int (mod), interp};
}
//------------------------------------------------------------------------------
// block delays have to be at least one block long. It has an additional
// parameter to overallocate memory to allow block processing "extra_spls".
//
// If the last element of a reverb loop is a delay of at least one block, all
// the outputs for a block can be fetched from existing previous work, before
// even starting processing the incoming sample batch.
//
// If those samples are to be used as feedback to sum with the current input
// batch, then the feedback samples need the last output from the past block
// to be added with the first incoming sample.
//
// So in total it is needed to fetch 1 block plus one old sample for feedback
// purposes, that's what the "extra_spls" parameter is for. This is used in
// conjunction with a value of 1 on "negative_offset" on the fetch_block
// function plus passing a buffer of the desired size + 1 spl.
static constexpr stage_data make_block_delay (u16 spls, u16 extra_spls = 1)
{
  return block_delay_data {
    fixpt_spls::from_int (spls), fixpt_spls::from_int (extra_spls)};
}
//------------------------------------------------------------------------------
template <class... Ts>
static constexpr stage_data make_multitap_delay (Ts&&... delay_spls)
{
  constexpr uint nargs = sizeof...(Ts);
  static_assert (nargs > 0);
  static_assert (nargs <= multitap_delay_data {}.spls.size());

  multitap_delay_data ret {};
  auto                tpl = std::forward_as_tuple (delay_spls...);
  mp11::mp_for_each<mp11::mp_iota_c<nargs>> ([&] (auto i) {
    ret.spls[i] = fixpt_spls::from_int (std::get<i> (tpl));
  });
  ret.count = nargs;
  return ret;
}
//------------------------------------------------------------------------------
template <class... Ts>
static constexpr stage_data make_parallel_delay (
  uint n_outs,
  Ts&&... spl_gain_pair)
{
  constexpr uint nargs = sizeof...(Ts);
  static_assert ((nargs % 2) == 0);
  static_assert (nargs > 0);
  static_assert ((nargs / 2) <= parallel_delay_data {}.spls.size());

  parallel_delay_data ret {};
  auto                tpl = std::forward_as_tuple (spl_gain_pair...);
  mp11::mp_for_each<mp11::mp_iota_c<nargs>> ([&] (auto i) {
    if constexpr ((i % 2) == 0) {
      ret.spls[i / 2] = fixpt_spls::from_int (std::get<i> (tpl));
    }
    else {
      ret.g[i / 2] = fixpt_t::from_float (std::get<i> (tpl));
    }
  });
  ret.count  = nargs / 2;
  ret.n_outs = n_outs;
  return ret;
}
//------------------------------------------------------------------------------
static constexpr stage_data make_hp (float g = 0)
{
  return filter_data {fixpt_t::from_float (g), false};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_lp (float g = 0)
{
  return filter_data {fixpt_t::from_float (g), true};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_crossover (
  float f    = 0,
  float g_lo = 0,
  float g_hi = 0)
{
  crossover_data r {};
  r.g[0]         = fixpt_t::from_float (f);
  r.band_gain[0] = fixpt_t::from_float (g_lo);
  r.band_gain[1] = fixpt_t::from_float (g_hi);
  r.n_bands      = 1;
  return r;
}
//------------------------------------------------------------------------------
static constexpr stage_data make_crossover2 (
  float f_lo  = 0,
  float f_hi  = 0,
  float g_lo  = 0,
  float g_mid = 0,
  float g_hi  = 0)
{
  crossover_data r {};
  r.g[0]         = fixpt_t::from_float (f_lo);
  r.g[1]         = fixpt_t::from_float (f_hi);
  r.band_gain[0] = fixpt_t::from_float (g_lo);
  r.band_gain[1] = fixpt_t::from_float (g_mid);
  r.band_gain[2] = fixpt_t::from_float (g_hi);
  r.n_bands      = 2;
  return r;
}
//------------------------------------------------------------------------------
static constexpr stage_data make_quantizer()
{
  return quantizer_data {};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_free_storage (uint n_elems)
{
  return free_storage_data {n_elems};
}
//------------------------------------------------------------------------------
// regular delays user thiran for small interpolation. This uses quantization.
// TODO: lerp
static constexpr stage_data make_variable_delay (
  uint min_spls,
  uint max_spls,
  bool abs_mod = false)
{
  return variable_delay_data {
    fixpt_spls::from_int (min_spls), fixpt_spls::from_int (max_spls), abs_mod};
}
//------------------------------------------------------------------------------
template <class Algorithm>
class spec_access {
public:
  //----------------------------------------------------------------------------
  static constexpr bool is_allpass (uint i)
  {
    return std::holds_alternative<allpass_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_comb (uint i)
  {
    return std::holds_alternative<comb_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_block_delay (uint i)
  {
    return std::holds_alternative<block_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_1tap_delay (uint i)
  {
    return std::holds_alternative<delay_data> (spec[i]) || is_block_delay (i);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_multitap_delay (uint i)
  {
    return std::holds_alternative<multitap_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_parallel_delay (uint i)
  {
    return std::holds_alternative<parallel_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_variable_delay (uint i)
  {
    return std::holds_alternative<variable_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool uses_absolute_delay_modulation (uint i)
  {
    if (is_variable_delay (i)) {
      return std::get<variable_delay_data> (spec[i]).mod_is_abs;
    }
    else {
      return false;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_filter (uint i)
  {
    return std::holds_alternative<filter_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_quantizer (uint i)
  {
    return std::holds_alternative<quantizer_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_lowpass (uint i)
  {
    if (is_filter (i)) {
      return std::get<filter_data> (spec[i]).is_lowpass;
    }
    return false;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_highpass (uint i)
  {
    if (is_filter (i)) {
      return !std::get<filter_data> (spec[i]).is_lowpass;
    }
    return false;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_crossover_n_bands (uint i)
  {
    if (std::holds_alternative<crossover_data> (spec[i])) {
      return std::get<crossover_data> (spec[i]).n_bands;
    }
    else {
      return 0;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_free_storage (uint i)
  {
    return std::holds_alternative<free_storage_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_free_storage_count (uint i)
  {
    if (is_free_storage (i)) {
      return std::get<free_storage_data> (spec[i]).count;
    }
    return 0;
  }
  //----------------------------------------------------------------------------
  static constexpr crossover_data get_crossover_data (uint i)
  {
    if (get_crossover_n_bands (i)) {
      return std::get<crossover_data> (spec[i]);
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr xspan<fixpt_t const> get_gains (uint i)
  {
    if (is_parallel_delay (i)) {
      auto& v = std::get<parallel_delay_data> (spec[i]);
      return {v.g.data(), v.count};
    }
    else if (is_allpass (i)) {
      return {&std::get<allpass_data> (spec[i]).g, 1};
    }
    else if (get_crossover_n_bands (i)) {
      return {get_crossover_data (i).g};
    }
    else if (is_filter (i)) {
      return {&std::get<filter_data> (spec[i]).g, 1};
    }
    else if (is_comb (i)) {
      return {&std::get<comb_data> (spec[i]).g, 1};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_t get_gain (uint i)
  {
    auto gain = get_gains (i);
    if (gain.size() == 1) {
      return gain[0];
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr xspan<fixpt_spls const> get_delays_spls (uint i)
  {
    if (is_multitap_delay (i)) {
      auto& v = std::get<multitap_delay_data> (spec[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_parallel_delay (i)) {
      auto& v = std::get<parallel_delay_data> (spec[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_allpass (i)) {
      return {&std::get<allpass_data> (spec[i]).spls, 1};
    }
    if (is_comb (i)) {
      return {&std::get<comb_data> (spec[i]).spls, 1};
    }
    else if (is_block_delay (i)) {
      return {&std::get<block_delay_data> (spec[i]).spls, 1};
    }
    else if (is_1tap_delay (i)) {
      return {&std::get<delay_data> (spec[i]).spls, 1};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_spls get_delay_spls (uint i)
  {
    auto spls = get_delays_spls (i);
    if (spls.size() == 1) {
      return spls[0];
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_spls_mod get_delay_mod_spls (uint i)
  {
    if (is_allpass (i)) {
      return std::get<allpass_data> (spec[i]).mod;
    }
    else if (is_comb (i)) {
      return std::get<comb_data> (spec[i]).mod;
    }
    else if (is_1tap_delay (i) && !is_block_delay (i)) {
      return std::get<delay_data> (spec[i]).mod;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_extra_spls (uint i)
  {
    if (is_block_delay (i)) {
      return std::get<block_delay_data> (spec[i]).extra_spls.to_int();
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr bool has_modulated_delay (uint i)
  {
    return get_delay_mod_spls (i).value() != 0;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_max_delay_spls (uint i)
  {
    if (is_parallel_delay (i) || is_multitap_delay (i)) {
      auto ds = get_delays_spls (i);
      return std::max_element (ds.begin(), ds.end())->to_int();
    }
    else if (is_variable_delay (i)) {
      return std::get<variable_delay_data> (spec[i]).max_spls.to_int();
    }
    else {
      return get_delay_spls (i).to_int() + get_delay_mod_spls (i).to_int()
        + ((has_modulated_delay (i)) ? 1 // extra tail spl for thiran/lerp
                                     : 0);
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_min_delay_spls (uint i)
  {
    if (is_parallel_delay (i) || is_multitap_delay (i)) {
      auto ds = get_delays_spls (i);
      return std::min_element (ds.begin(), ds.end())->to_int();
    }
    if (is_variable_delay (i)) {
      return std::get<variable_delay_data> (spec[i]).min_spls.to_int();
    }
    else {
      return get_delay_spls (i).to_int() - get_delay_mod_spls (i).to_int();
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_n_outs (uint i)
  {
    if (is_parallel_delay (i)) {
      return std::get<parallel_delay_data> (spec[i]).n_outs;
    }
    else {
      return 1;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr std::optional<interpolation> get_interp (uint i)
  {
    if (is_allpass (i) && has_modulated_delay (i)) {
      return std::get<allpass_data> (spec[i]).interp;
    }
    else if (is_comb (i) && has_modulated_delay (i)) {
      return std::get<comb_data> (spec[i]).interp;
    }
    else if (is_1tap_delay (i) && has_modulated_delay (i)) {
      return std::get<delay_data> (spec[i]).interp;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr std::size_t size() { return spec.size(); }
  //----------------------------------------------------------------------------
private:
  static constexpr auto spec {Algorithm::get_spec()};
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
namespace state {
// these need to be outside of engine, as the constexpr function has to be
// fully defined before.
//------------------------------------------------------------------------------
struct empty {
  // struct clang_bug_workaround {};
  // clang_bug_workaround v;
};
//------------------------------------------------------------------------------
struct y1_float {
  float y1;
};
//------------------------------------------------------------------------------
template <uint N>
struct y1_float_arr {
  std::array<y1_float, N> arr;
};
//------------------------------------------------------------------------------
struct y1_fixpt {
  fixpt_t y1;
  s16     y1_err;
};
//------------------------------------------------------------------------------
template <uint N>
struct y1_fixpt_arr {
  std::array<y1_fixpt, N> arr;
};

//------------------------------------------------------------------------------
struct allpass_fixpt {
  s16 x_err;
  s16 u_err;
};
//------------------------------------------------------------------------------
struct quantizer_fixpt {
  s16 err;
};
//------------------------------------------------------------------------------
template <uint N>
struct quantizer_arr_fixpt {
  std::array<s16, N> err;
};
//------------------------------------------------------------------------------
struct allpass_and_y1_fixpt : public allpass_fixpt, public y1_fixpt {};
//------------------------------------------------------------------------------
struct quantizer_and_y1_fixpt : public quantizer_fixpt, public y1_fixpt {};
//------------------------------------------------------------------------------
template <class value_type, uint N>
struct free_storage {
  std::array<value_type, N> sto;
};
//------------------------------------------------------------------------------
template <class value_type, uint Idx, class SpecAccess>
static constexpr auto get_state_type()
{
  constexpr auto interp = SpecAccess::get_interp (Idx);
  if constexpr (is_fixpt_v<value_type>) {
    if constexpr (SpecAccess::is_allpass (Idx)) {
      if constexpr (!!interp && *interp == interpolation::thiran) {
        return allpass_and_y1_fixpt {};
      }
      else {
        return allpass_fixpt {};
      }
    }
    else if constexpr (SpecAccess::is_comb (Idx)) {
      if constexpr (!!interp && *interp == interpolation::thiran) {
        return quantizer_and_y1_fixpt {};
      }
      else {
        return quantizer_fixpt {};
      }
    }
    else if constexpr (SpecAccess::is_1tap_delay (Idx)) {
      if constexpr (!!interp && *interp == interpolation::thiran) {
        return y1_fixpt {};
      }
      else {
        return empty {};
      }
    }
    else if constexpr (SpecAccess::is_filter (Idx)) {
      return y1_fixpt {};
    }
    else if constexpr (SpecAccess::get_crossover_n_bands (Idx)) {
      return y1_fixpt_arr<SpecAccess::get_crossover_n_bands (Idx)> {};
    }
    else if constexpr (SpecAccess::is_quantizer (Idx)) {
      return quantizer_fixpt {};
    }
    else if constexpr (SpecAccess::is_parallel_delay (Idx)) {
      return quantizer_arr_fixpt<SpecAccess::get_n_outs (Idx)> {};
    }
    else if constexpr (SpecAccess::is_free_storage (Idx)) {
      return free_storage<
        value_type,
        SpecAccess::get_free_storage_count (Idx)> {};
    }
    else {
      return empty {};
    }
  }
  else {
    static_assert (std::is_same_v<value_type, float>);
    if constexpr (
      (!!interp && *interp == interpolation::thiran)
      || SpecAccess::is_filter (Idx)) {
      return y1_float {};
    }
    else if constexpr (SpecAccess::get_crossover_n_bands (Idx)) {
      return y1_float_arr<SpecAccess::get_crossover_n_bands (Idx)> {};
    }
    else if constexpr (SpecAccess::is_free_storage (Idx)) {
      return free_storage<
        value_type,
        SpecAccess::get_free_storage_count (Idx)> {};
    }
    else {
      return empty {};
    }
  }
}
//------------------------------------------------------------------------------
template <class value_type, class SpecAccess>
struct index_to_state_qfn {

  template <class Idx>
  using fn = decltype (get_state_type<value_type, Idx::value, SpecAccess>());
};
//------------------------------------------------------------------------------

} // namespace state

namespace delay {

//------------------------------------------------------------------------------
// TODO: move to delay_line? this one provides the functions to get the segments
template <class T, uint Size>
class static_buffer {
public:
  //----------------------------------------------------------------------------
  using value_type           = T;
  static constexpr uint size = Size;
  //----------------------------------------------------------------------------
  void insert (value_type v)
  {
    // TODO: should the advance should be done before? to get this sample now
    // requires a call to read spl0
    _z[_pos] = v;
    advance_one();
  }
  //----------------------------------------------------------------------------
  // returns the memory for a block a insertion, updates current position
  std::array<xspan<value_type>, 2> insert_block (uint blocksize)
  {
    assert (blocksize);
    assert (_z.size() >= blocksize);

    uint block1 = _pos;
    uint end    = block1 + blocksize - 1;
    end -= (end >= _z.size()) ? _z.size() : 0;
    _pos = end;
    advance_one();

    uint block1sz, block2sz;
    if (block1 <= end) {
      // contiguous, no wraparound
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = _z.size() - block1;
      block2sz = blocksize - block1sz;
    }
    return {xspan {&_z[block1], block1sz}, xspan {&_z[0], block2sz}};
  }
  //----------------------------------------------------------------------------
  value_type read (uint delay_spls)
  {
    assert (delay_spls <= _z.size());
    uint z = _pos - delay_spls;
    z += (z >= _z.size()) ? _z.size() : 0;
    return _z[z];
  }
  //----------------------------------------------------------------------------
  std::array<xspan<value_type>, 2> read_block (uint blocksize, uint del_spls)
  {
    assert (blocksize);
    assert (del_spls <= _z.size());
    assert (del_spls >= blocksize);
    assert (_z.size() >= blocksize);

    uint block1 = _pos - del_spls;
    block1 += (block1 >= _z.size()) ? _z.size() : 0;
    uint end = _pos - del_spls + blocksize - 1;
    end += (end >= _z.size()) ? _z.size() : 0;

    uint block1sz, block2sz;
    if (block1 <= end) {
      // contiguous
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = _z.size() - block1;
      block2sz = blocksize - block1sz;
    }
    return {xspan {&_z[block1], block1sz}, xspan {&_z[0], block2sz}};
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void advance_one()
  {
    ++_pos;
    if (_pos == _z.size()) {
      _pos = 0;
    }
  }
  //----------------------------------------------------------------------------
  std::array<value_type, size> _z {};
  uint                         _pos {};
};
//------------------------------------------------------------------------------

template <class T, uint Size, class C>
class static_encoded_buffer {
public:
  //----------------------------------------------------------------------------
  using value_type     = T;
  using enconder       = C;
  using enconding_type = typename enconder::value_type;
  static constexpr bool has_enconding
    = !std::is_same_v<value_type, enconding_type>;
  static constexpr uint size = Size;
  //----------------------------------------------------------------------------
  value_type read (uint spl)
  {
    assert (spl);
    if constexpr (has_enconding) {
      return enconder::decode (_z.read (spl));
    }
    else {
      return _z.read (spl);
    }
  }
  //----------------------------------------------------------------------------
  void push (value_type v)
  {
    if constexpr (has_enconding) {
      _z.insert (enconder::encode (v));
    }
    else {
      _z.insert (v);
    }
  }
  //----------------------------------------------------------------------------
  void read_block (xspan<value_type> dst, uint del_spls)
  {
    std::array<xspan<enconding_type>, 2> src
      = _z.read_block (dst.size(), del_spls);
    assert (dst.size() == (src[0].size() + src[1].size()));
    if constexpr (has_enconding) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < src[0].size(); ++i) {
        dst[i] = enconder::decode (src[0][i]);
      }
      dst.cut_head (src[0].size());
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < src[1].size(); ++i) {
        dst[i] = enconder::decode (src[1][i]);
      }
    }
    else {
      xspan_memdump (dst.data(), src[0]);
      dst.cut_head (src[0].size());
      xspan_memdump (dst.data(), src[1]);
    }
  }
  //----------------------------------------------------------------------------
  void push_block (xspan<value_type const> src)
  {
    std::array<xspan<enconding_type>, 2> dst = _z.insert_block (src.size());
    assert (src.size() == (dst[0].size() + dst[1].size()));
    if constexpr (has_enconding) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst[0].size(); ++i) {
        dst[0][i] = enconder::encode (src[i]);
      }
      src.cut_head (dst[0].size());
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst[1].size(); ++i) {
        dst[1][i] = enconder::encode (src[i]);
      }
    }
    else {
      xspan_memcpy (dst[0], src);
      src.cut_head (dst[0].size());
      xspan_memcpy (dst[1], src);
    }
  }
  //----------------------------------------------------------------------------
private:
  static_buffer<enconding_type, size> _z {};
};
//------------------------------------------------------------------------------

template <class T>
struct no_enconding {
  using value_type = T;
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
enum class data_type : uint { fixpt16, float16, float32 };
//------------------------------------------------------------------------------
template <data_type Te, uint Size>
class delay_buffer;

template <>
class delay_buffer<data_type::fixpt16, 0> {};

template <>
class delay_buffer<data_type::float16, 0> {};

template <>
class delay_buffer<data_type::float32, 0> {};

template <uint Size>
class delay_buffer<data_type::fixpt16, Size>
  : public static_encoded_buffer<fixpt_t, Size, no_enconding<fixpt_t>> {
public:
  using value_type = fixpt_t;
};

using float16_encoder = f16pack<5, -1, f16pack_dftz | f16pack_clamp>;

template <uint Size>
class delay_buffer<data_type::float16, Size>
  : public static_encoded_buffer<float, Size, float16_encoder> {
public:
  using value_type =
    typename static_encoded_buffer<float, Size, float16_encoder>::value_type;
};

template <uint Size>
class delay_buffer<data_type::float32, Size>
  : public static_encoded_buffer<float, Size, no_enconding<float>> {
public:
  using value_type = float;
};
//------------------------------------------------------------------------------
template <data_type Te, class SpecAccess>
struct index_to_delay_buffer_qfn {

  template <class Idx>
  using fn = delay_buffer<
    Te,
    (uint) SpecAccess::get_max_delay_spls (Idx::value)
      + (uint) SpecAccess::get_delay_extra_spls (Idx::value)>;
};
//------------------------------------------------------------------------------
} // namespace delay

struct lofiverb_io_tag {};

struct defaulted_tag {};
static constexpr defaulted_tag defaulted {};

struct add_to_out_tag : public lofiverb_io_tag {};
static constexpr add_to_out_tag add_to {};

struct overwrite_out_tag : public lofiverb_io_tag {};
static constexpr overwrite_out_tag overwrite {};

//----------------------------------------------------------------------------
template <uint... Idxs>
struct stage_list {}; // stage list, used for avoiding the template keyword
//------------------------------------------------------------------------------
// A class to abstract 16-bit storage, queue access and common DSP operations
// when building reverbs based on allpass loops. Both on fixed and floating
// point.
template <class Algorithm, delay::data_type Data_type, uint Max_block_size>
class engine {
public:
  //----------------------------------------------------------------------------
  using value_type = typename delay::delay_buffer<Data_type, 1>::value_type;
  using algorithm  = Algorithm;
  using spec       = spec_access<Algorithm>;
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = Max_block_size;
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes() { return sizeof *_delay; }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<u8> mem) // needs malloc/new alignment
  {
    assert (mem.size() >= sizeof (*_delay));
    // Tuple of trivial types, no destructor required.
    _delay = new (mem.data()) std::decay_t<decltype (*_delay)> {};
  }
  //----------------------------------------------------------------------------
  template <uint Idx, uint... Idxs>
  auto get_gain_for_rt60 (stage_list<Idx, Idxs...>, float t_sec, uint srate)
  {
    constexpr uint vec_size    = 4;
    using V                    = vec<float, vec_size>;
    constexpr uint  n_vals     = 1 + sizeof...(Idxs);
    constexpr uint  n_vecs     = div_ceil (n_vals, vec_size);
    constexpr uint  n_rounded  = n_vecs * vec_size;
    constexpr float log2_milli = -9.965784284662087f; // log2(0.001)

    alignas (V) std::array<float, n_rounded> spls {
      {(float) spec::get_delay_spls (Idx),
       (float) spec::get_delay_spls (Idxs)...}};

    std::array<float, n_vals> ret_flt;

    for (uint i = 0; i < n_vecs; ++i) {
      auto exps = vec_load<V> (&spls[i * vec_size]);
      exps /= srate * t_sec;
      auto gains = vec_exp2_2dat (exps * vec_set<V> (log2_milli));
      if (i != (n_vecs - 1) || (n_rounded == n_vals)) {
        vec_store_unaligned (&ret_flt[i * vec_size], gains);
      }
      else {
        for (uint j = 0; j < (n_vals % vec_size); ++j) {
          ret_flt[i * vec_size + j] = gains[j];
        }
      }
    }
    if constexpr (std::is_same_v<value_type, float>) {
      return ret_flt;
    }
    else if constexpr (is_fixpt_v<value_type>) {
      return std::apply (
        [] (auto... x) { return make_array (value_type::from_float (x)...); },
        ret_flt);
    }
    else {
      static_assert (Idx != Idx, "unsupported type");
    }
  }
  //----------------------------------------------------------------------------
  // Pure delays with "size > blocksize" are placed before feedback loops to
  // enable block processing, so it is possible to fetch the future
  // feedbacks/outputs (minus the first) at once and then to push them all at
  // once. This is exactly what fetch/push accomplish.
  //
  // dst has to contain one element more at the head, where the previous
  // output will be placed. Once the feedback is applied to a current input,
  // this head feedback sample (last output) can be dropped.
  template <uint Idx, std::enable_if_t<spec::is_block_delay (Idx)>* = nullptr>
  void fetch_block (
    stage_list<Idx>,
    xspan<value_type> dst,
    uint              negative_offset = 0)
  {
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    std::get<Idx> (*_delay).read_block (
      dst, spec::get_delay_spls (Idx).to_int() + negative_offset);
  }
  //----------------------------------------------------------------------------
  // Fetch samples from the queue of a processed allpass. The allpass has to be
  // independenly processed.
  template <uint Idx, std::enable_if_t<spec::is_allpass (Idx)>* = nullptr>
  void fetch_block (stage_list<Idx>, xspan<value_type> dst, uint spls)
  {
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (dst);
    assert (spls <= spec::get_max_delay_spls (Idx));
    assert (spls >= max_block_size);
    std::get<Idx> (*_delay).read_block (dst, spls);
  }
  //----------------------------------------------------------------------------
  // fetches the output and feedback from an allpass. The output will
  // be able to be used directly. The feedback signal will have to be manually
  // pushed (see push).
  //
  // The intent is to be able to process further (no gain) the feedback
  // signal, e.g. filter it, before insertion.
  //
  // io[in, out = input samples, outoput samples at out
  // fb[out] signal, the feedback. Can be processed (e.g) filtered and the a
  // call to push is required
  template <uint Idx, class Lfo, class G>
  void fetch_block (
    stage_list<Idx>,
    xspan<value_type>       out,
    xspan<value_type>       fb,
    xspan<value_type const> in,
    Lfo&&                   lfo,
    G&&                     gain)
  {
    static_assert (spec::is_allpass (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (out && fb && in);
    assert (fb.size() >= in.size());
    assert (out.size() >= in.size());

    auto g_gen = get_gain_generator<Idx> (std::forward<G> (gain));

    // no overlap, can run block-wise
    fetch_block_optmod<Idx> (
      xspan {fb.data(), in.size()}, std::forward<Lfo> (lfo));
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      if constexpr (spec::is_allpass (Idx)) {
        std::tie (out[i], fb[i]) = run_allpass<Idx> (in[i], fb[i], g_gen (i));
      }
      else {
        std::tie (out[i], fb[i]) = run_comb<Idx> (in[i], fb[i], g_gen (i));
      }
    }
  }

  template <uint Idx, class Lfo, class G>
  void fetch_block (
    stage_list<Idx>,
    xspan<value_type> io,
    xspan<value_type> fb,
    Lfo&&             lfo,
    G&&               gain)
  {
    static_assert (spec::is_allpass (Idx));
    fetch_block<Idx> (
      io, fb, io.to_const(), std::forward<Lfo> (lfo), std::forward<G> (gain));
  }
  //----------------------------------------------------------------------------
  // Fetch for combs, gets the delay modulated and with gain applied
  template <uint Idx, class value_type, class Lfo, class G>
  void fetch_block (stage_list<Idx>, xspan<value_type> fb, Lfo&& lfo, G&& gain)
  {
    static_assert (spec::is_comb (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (fb);
    fetch_block_optmod<Idx> (fb, std::forward<Lfo> (lfo));
    auto g_gen = get_gain_generator<Idx> (std::forward<G> (gain));
    run_quantizer<Idx> (fb.data(), fb, [&] (auto v, uint i) {
      return v * g_gen (i);
    });
  }
  // see comment on fetch_block overloads
  //----------------------------------------------------------------------------
  // Push for combs, gets the feedback signal (obtained from fetch_block and
  // maybe filtered) and the input, returns the sum.
  template <uint Idx>
  void push (
    stage_list<Idx>,
    xspan<value_type>       sum,
    xspan<value_type const> fb,
    xspan<value_type const> in)
  {
    static_assert (spec::is_comb (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (sum && fb && in);
    assert (fb.size() && in.size());
    span_add (sum, fb, in);
    std::get<Idx> (*_delay).push_block (xspan {sum.data(), in.size()});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class value_type>
  void push (stage_list<Idx>, xspan<value_type const> src)
  {
    static_assert (spec::is_block_delay (Idx) || spec::is_allpass (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    std::get<Idx> (*_delay).push_block (src);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  xspan<value_type> get_storage (stage_list<Idx>)
  {
    static_assert (spec::is_free_storage (Idx));
    return {std::get<Idx> (_states).sto};
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void run (
    stage_list<Idx>,
    xspan<value_type>       out,
    xspan<value_type const> in,
    Ts&&... args)
  {
    assert (out && in);
    assert (out.size() >= in.size());
    if constexpr (
      spec::is_allpass (Idx) || spec::is_comb (Idx)
      || spec::is_1tap_delay (Idx)) {
      constexpr uint n_args          = sizeof...(Ts);
      constexpr uint expected_n_args = 2;

      if constexpr (expected_n_args == n_args) {
        run_ap_comb_or_delay<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
      else {
        // add null on the remaining positions
        using filler_defaulted = mp11::
          mp_repeat_c<std::tuple<defaulted_tag>, expected_n_args - n_args>;
        std::apply (
          [&] (auto&&... targs) {
            run_ap_comb_or_delay<Idx> (
              out.data(), in, std::forward<decltype (targs)> (targs)...);
          },
          std::tuple_cat (
            std::forward_as_tuple (std::forward<Ts> (args)...),
            filler_defaulted {}));
      }
    }
    else if constexpr (spec::is_variable_delay (Idx)) {
      static_assert (
        sizeof...(args) == 1, "variable delays always require mod");
      run_variable_delay<Idx> (out, in, std::forward<Ts> (args)...);
    }
    else if constexpr (spec::is_filter (Idx)) {
      if constexpr (spec::is_lowpass (Idx)) {
        run_lp<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
      else {
        run_hp<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
    }
    else if constexpr (spec::get_crossover_n_bands (Idx)) {
      run_crossover<Idx> (out.data(), in, std::forward<Ts> (args)...);
    }
    else if constexpr (spec::is_quantizer (Idx)) {
      run_quantizer<Idx> (out.data(), in, std::forward<Ts> (args)...);
    }
    else {
      static_assert (Idx != Idx, "Invalid");
    }
  }
  //----------------------------------------------------------------------------
  template <
    uint Idx,
    class U,
    class Tag,
    class... Ts,
    std::enable_if_t<std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<Tag>>>* = nullptr>
  void run (
    stage_list<Idx>,
    xspan<value_type const> in,
    Tag                     t,
    U&&                     out,
    Ts&&... outs)
  {
    assert (in);
    if constexpr (spec::is_multitap_delay (Idx)) {
      run_multitap_delay<Idx> (
        in, t, std::forward<U> (out), std::forward<Ts> (outs)...);
    }
    else if constexpr (spec::is_parallel_delay (Idx)) {
      run_parallel_delays<Idx> (
        in, t, std::forward<U> (out), std::forward<Ts> (outs)...);
    }
  }

  template <
    uint Idx,
    class Tag,
    std::enable_if_t<std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<Tag>>>* = nullptr>
  void run (stage_list<Idx> l, xspan<value_type> io, Tag t)
  {
    run (l, io.to_const(), t, io);
  }

  template <
    uint Idx,
    class U,
    class... Ts,
    std::enable_if_t<!std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<U>>>* = nullptr>
  void run (stage_list<Idx> l, xspan<value_type> io, U&& arg1, Ts&&... args)
  {
    run (
      l, io, io.to_const(), std::forward<U> (arg1), std::forward<Ts> (args)...);
  }

  template <uint Idx>
  void run (stage_list<Idx> l, xspan<value_type> io)
  {
    run (l, io, io.to_const());
  }
  //----------------------------------------------------------------------------
  // for arbitrarily nested allpasses or crossovers
  template <uint Idx1, uint Idx2, uint... Idxs, class... Ts>
  void run (
    stage_list<Idx1, Idx2, Idxs...>,
    xspan<value_type>       out,
    xspan<value_type const> in,
    Ts&&... args)
  {
    if constexpr (
      spec::is_allpass (Idx1) && can_be_placed_on_nested_allpass<Idx2>
      && (... && can_be_placed_on_nested_allpass<Idxs>) ) {
      run_nested_ap<Idx1, Idx2, Idxs...> (
        out.data(), in, std::forward<Ts> (args)...);
    }
    else {
      static_assert (Idx1 != Idx1, "Invalid type on one of the indexes");
    }
  }

  template <uint Idx1, uint Idx2, uint... Idxs, class... Ts>
  void run (
    stage_list<Idx1, Idx2, Idxs...> l,
    xspan<value_type>               io,
    Ts&&... args)
  {
    run (l, io, io.to_const(), std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
private:
  template <uint Idx>
  static constexpr bool can_be_placed_on_nested_allpass = spec::is_allpass (Idx)
    || spec::is_1tap_delay (Idx) || spec::is_lowpass (Idx)
    || spec::is_highpass (Idx) || spec::get_crossover_n_bands (Idx) != 0;
  //----------------------------------------------------------------------------
  template <class Lfo>
  static constexpr auto get_generic_generator (Lfo&& v)
  {
    using Lfo_no_cv_ref = std::remove_cv_t<std::remove_reference_t<Lfo>>;
    if constexpr (is_generator<Lfo_no_cv_ref>()) {
      return std::forward<Lfo_no_cv_ref> (v);
    }
    else if constexpr (is_array_subscriptable_v<Lfo_no_cv_ref>) {
      return [v] (uint i) { return v[i]; };
    }
    else {
      return defaulted_tag {};
    }
  }
  //----------------------------------------------------------------------------
  template <class Lfo>
  static constexpr auto get_sample_generator (Lfo&& v)
  {
    using Lfo_no_cv_ref = std::remove_cv_t<std::remove_reference_t<Lfo>>;
    if constexpr (is_generator<Lfo_no_cv_ref>()) {
      return std::forward<Lfo_no_cv_ref> (v);
    }
    else if constexpr (is_uint_generator<Lfo_no_cv_ref>()) {
      return std::forward<Lfo_no_cv_ref> (v);
    }
    else if constexpr (is_array_subscriptable_v<Lfo_no_cv_ref>) {
      return [v] (uint i) { return v[i]; };
    }
    else {
      return defaulted_tag {};
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Gain>
  static constexpr auto get_gain_generator (Gain&& v)
  {
    using Gain_no_cv_ref = std::remove_cv_t<std::remove_reference_t<Gain>>;
    if constexpr (is_generator<Gain_no_cv_ref>()) {
      return std::forward<Gain_no_cv_ref> (v);
    }
    else if constexpr (is_array_subscriptable_v<Gain_no_cv_ref>) {
      return [v] (uint i) { return v[i]; };
    }
    else if constexpr (std::is_same_v<value_type, Gain_no_cv_ref>) {
      return [v] (uint) { return v; };
    }
    else if constexpr (std::is_same_v<defaulted_tag, Gain_no_cv_ref>) {
      return [] (uint) { return (value_type) spec::get_gain (Idx); };
    }
    else {
      static_assert (
        sizeof (Gain_no_cv_ref) != sizeof (Gain_no_cv_ref), "unknown type");
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  value_type get (uint delay_spls)
  {
    return std::get<Idx> (*_delay).read (delay_spls);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  value_type read_next()
  {
    constexpr auto spls = spec::get_delay_spls (Idx).to_int();
    return get<Idx> (spls);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Func>
  void run_quantizer (value_type* out, xspan<value_type const> in, Func&& fn)
  {
    if constexpr (is_fixpt_v<value_type>) {
      // decay with error-feedback/ fraction saving
      state::quantizer_fixpt& st  = std::get<Idx> (_states);
      auto                    err = st.err;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        auto [v, err_] = round (fn (in[i], i), err);
        out[i]         = v;
        err            = err_;
      }
      st.err = err;
    }
    else {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        out[i] = fn (in[i], i);
      }
    }
  }
  //----------------------------------------------------------------------------
  enum onepole_type { lp, hp };
  //----------------------------------------------------------------------------
  template <onepole_type Type, class Gain, class State>
  void run_1pole (
    value_type*             out,
    xspan<value_type const> in,
    Gain                    g,
    State&                  st)
  {
    assert (out && in);

    if constexpr (is_fixpt_v<value_type>) {
      auto y1  = st.y1;
      auto err = st.y1_err;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        auto         x = (fixpt_acum_t) in[i];
        fixpt_acum_t gw;
        if constexpr (is_generator<Gain>()) {
          gw = (fixpt_acum_t) g (i);
        }
        else if constexpr (is_array_subscriptable_v<Gain>) {
          gw = (fixpt_acum_t) g[i];
        }
        else {
          gw = (fixpt_acum_t) g;
        }
        auto y           = (1_r - g) * x + g * y1;
        auto [y1_, err_] = truncate (y, err);
        y1               = y1_;
        err              = err_;
        if constexpr (Type == onepole_type::lp) {
          out[i] = y1;
        }
        else {
          static_assert (Type == onepole_type::hp);
          out[i] = (fixpt_t) (x - y1);
        }
      }
      st.y1     = y1;
      st.y1_err = err;
    }
    else {
      auto y1 = st.y1;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        auto       x = in[i];
        value_type gv;
        if constexpr (is_generator<Gain>()) {
          gv = (value_type) g (i);
        }
        else if constexpr (is_array_subscriptable_v<Gain>) {
          gv = (value_type) g[i];
        }
        else {
          gv = (value_type) g;
        }
        y1 = (1.f - gv) * x + gv * y1;
        if constexpr (Type == onepole_type::lp) {
          out[i] = y1;
        }
        else {
          static_assert (Type == onepole_type::hp);
          out[i] = (value_type) (x - y1);
        }
      }
      st.y1 = y1;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_lp (value_type* out, xspan<value_type const> in, U&& g)
  {
    run_1pole<onepole_type::lp> (
      out, in, std::forward<U> (g), std::get<Idx> (_states));
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void run_lp (value_type* out, xspan<value_type const> in)
  {
    run_lp<Idx> (out, in, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_hp (value_type* out, xspan<value_type const> in, U&& g)
  {
    run_1pole<onepole_type::hp> (
      out, in, std::forward<U> (g), std::get<Idx> (_states));
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void run_hp (value_type* out, xspan<value_type const> in)
  {
    run_hp<Idx> (out, in, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, std::size_t N>
  void run_crossover_impl (
    value_type*                   out,
    xspan<value_type const>       in,
    std::array<value_type, N>     freq,
    std::array<value_type, N + 1> gain)
  {
    constexpr uint                   n_bands = N;
    std::array<block_array, n_bands> band_mem;

    auto& st = std::get<Idx> (_states);

    xspan_memdump (band_mem[0].data(), in);
    std::array<value_type*, n_bands + 1> band;
    for (uint i = 0; i < n_bands; ++i) {
      band[i] = band_mem[i].data();
    }
    band[n_bands] = out;
    // compute bands, scale each one except the lowpass
    for (uint b = 0; b < n_bands; ++b) {
      // starting from higher to lower freq.
      uint idx = n_bands - 1 - b;
      run_1pole<onepole_type::lp> (
        band[b + 1], xspan {band[b], in.size()}, freq[idx], st.arr[b]);

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        // finishing the previous band by subtracting this one and gain scaling
        band[b][i] = (value_type) (band[b][i] - band[b + 1][i]);
        band[b][i] = (value_type) (band[b][i] * gain[idx + 1]);
      }
    }
    // scale lowpass band and sum
    std::apply (
      [out, &band, &gain, &in] (auto&&... bands) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          band[n_bands][i] = (value_type) (band[n_bands][i] * gain[0]);
          out[i]           = (value_type) (bands[i] + ...);
        }
      },
      band);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_crossover (value_type* out, xspan<value_type const> in)
  {
    static constexpr auto                 d = spec::get_crossover_data (Idx);
    std::array<value_type, d.n_bands>     g;
    std::array<value_type, d.n_bands + 1> bg;
    for (uint i = 0; i < d.n_bands; ++i) {
      g[i]  = d.g[i];
      bg[i] = d.bg[i];
    }
    bg[d.n_bands] = d.bg[d.n_bands];
    run_crossover_impl<Idx> (out, in, g, bg);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_crossover (
    value_type*             out,
    xspan<value_type const> in,
    U&&                     f_lo,
    U&&                     g_lo,
    U&&                     g_hi)
  {
    run_crossover_impl<Idx> (
      out,
      in,
      make_array ((value_type) f_lo),
      make_array ((value_type) g_lo, (value_type) g_hi));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_crossover (
    value_type*             out,
    xspan<value_type const> in,
    U                       f_lo,
    U                       g_lo,
    U                       f_mid,
    U                       g_mid,
    U                       g_hi)
  {
    run_crossover_impl<Idx> (
      out,
      in,
      make_array ((value_type) f_lo, (value_type) f_mid),
      make_array ((value_type) g_lo, (value_type) g_mid, (value_type) g_hi));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_thiran (xspan<value_type> dst, LF&& lfo_gen)
  {
    static_assert (spec::has_modulated_delay (Idx));
    constexpr auto delay_spls     = spec::get_delay_spls (Idx).add_sign();
    constexpr auto delay_mod_spls = spec::get_delay_mod_spls (Idx).add_sign();
    auto&          st             = std::get<Idx> (_states);

    if constexpr (is_fixpt_v<value_type>) {
      using fixpt_th = decltype (fixpt_acum_t {}.resize<2, -2>());

      auto y1  = st.y1;
      auto err = st.y1_err;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        auto fixpt_spls = (delay_spls + (lfo_gen (i) * delay_mod_spls));
        auto n_spls     = fixpt_spls.to_int();
        n_spls -= i;
        auto n_spls_frac = (fixpt_th) fixpt_spls.fractional();
        auto d           = n_spls_frac + 0.418_r; // this might exceed 1
        auto a           = (1_r - d) / (1_r + d); // 0.4104 to -1
        auto z0          = (fixpt_th) get<Idx> (n_spls - 1);
        auto z1          = (fixpt_th) get<Idx> (n_spls);
        auto y           = (z0 * a) + z1 - (a * y1);
        auto [y1_, err_] = truncate (y, err);
        y1               = y1_;
        err              = err_;
        dst[i]           = y1;
      }
      st.y1     = y1;
      st.y1_err = err;
    }
    else {
      value_type y1 = st.y1;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        value_type fixpt_spls
          = delay_spls.to_floatp() + (lfo_gen (i) * delay_mod_spls.to_floatp());
        uint       n_spls      = (uint) fixpt_spls;
        value_type n_spls_frac = fixpt_spls - n_spls;
        n_spls -= i;

        value_type d = n_spls_frac + 0.418_r;
        value_type a = (1_r - d) / (1_r + d); // 0.4104 to -1

        auto z0 = get<Idx> (n_spls - 1);
        auto z1 = get<Idx> (n_spls);
        y1      = z0 * a + z1 - a * y1;
        dst[i]  = y1;
      }
      st.y1 = y1;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_lerp (xspan<value_type> dst, LF&& lfo_gen)
  {
    static_assert (spec::has_modulated_delay (Idx));
    constexpr auto delay_spls     = spec::get_delay_spls (Idx).add_sign();
    constexpr auto delay_mod_spls = spec::get_delay_mod_spls (Idx).add_sign();

    if constexpr (is_fixpt_v<value_type>) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        auto fixpt_spls = (delay_spls + (lfo_gen (i) * delay_mod_spls));
        auto n_spls     = fixpt_spls.to_int();
        n_spls -= i;
        auto n_spls_frac = fixpt_spls.fractional();
        auto zb          = get<Idx> (n_spls + 1);
        auto z           = get<Idx> (n_spls);
        static_assert (n_spls_frac.n_int == 0);
        dst[i]
          = (fixpt_t) (z * (n_spls_frac.max() - n_spls_frac) + zb * n_spls_frac);
      }
    }
    else {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        value_type fixpt_spls
          = delay_spls.to_floatp() + (lfo_gen (i) * delay_mod_spls.to_floatp());
        uint  n_spls      = (uint) fixpt_spls;
        float n_spls_frac = fixpt_spls - n_spls;
        n_spls -= i;
        auto zb = get<Idx> (n_spls + 1);
        auto z  = get<Idx> (n_spls);
        dst[i]  = z * (1.f - n_spls_frac) + zb * n_spls_frac;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  auto run_allpass (value_type in, value_type yn, U g)
  {
    static_assert (spec::is_allpass (Idx));
    auto u = in + yn * g;
    auto x = yn - u * g;
    if constexpr (std::is_floating_point_v<value_type>) {
      return std::make_tuple (x, u);
    }
    else {
      auto& st          = std::get<Idx> (_states);
      auto [q_x, x_err] = truncate (x, st.x_err);
      auto [q_u, u_err] = truncate (u, st.u_err);
      st.x_err          = x_err;
      st.u_err          = u_err;
      return std::make_tuple (q_x, q_u);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  auto run_comb (value_type in, value_type z, U g)
  {
    static_assert (spec::is_comb (Idx));
    auto u = in + z * g;
    if constexpr (std::is_floating_point_v<value_type>) {
      return std::make_tuple (z, u);
    }
    else {
      auto& st        = std::get<Idx> (_states.fix);
      auto [q_u, err] = truncate (u, st.err);
      st.err          = err;
      return std::make_tuple (z, q_u);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Lfo>
  void fetch_block_optmod (xspan<value_type> dst, Lfo&& lfo)
  {
    constexpr auto minsz = spec::get_min_delay_spls (Idx);
    assert (dst.size() <= minsz);
    if constexpr (minsz >= max_block_size) {
      constexpr std::optional<interpolation> interp = spec::get_interp (Idx);
      if constexpr (!interp) {
        std::get<Idx> (*_delay).read_block (
          dst, (uint) spec::get_delay_spls (Idx));
      }
      else {
        auto lfo_gen = get_generic_generator (std::forward<Lfo> (lfo));
        if constexpr (*interp == interpolation::linear) {
          run_lerp<Idx> (dst, lfo_gen);
        }
        else {
          static_assert (*interp == interpolation::thiran);
          run_thiran<Idx> (dst, lfo_gen);
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  // Template function to run:
  // - delays
  // - modulated delays
  // - non-nested allpasses
  // - modulated non-nested allpasses.
  //
  // Lfo and gain parameters can be either value generators (lambdas with
  // value_type(int) signature), spans or dummy types. Use "defaulted_tag" to
  // pass defaulted parameters (no delay modulation, same allpass gain as in the
  // specification)
  template <uint Idx, class Lfo, class G>
  void run_ap_comb_or_delay (
    value_type*             out,
    xspan<value_type const> in,
    Lfo&&                   lfo,
    G&&                     gain)
  {
    // reminder, out and in can alias
    xspan          sout {out, in.size()};
    constexpr auto minsz = spec::get_min_delay_spls (Idx);
    if constexpr (minsz >= max_block_size) {
      // no overlap, can run block-wise
      block_array z_mem;
      xspan       z {z_mem.data(), in.size()};
      if constexpr (spec::is_allpass (Idx) || spec::is_comb (Idx)) {
        fetch_block (
          stage_list<Idx> {},
          sout,
          z,
          in.to_const(),
          std::forward<Lfo> (lfo),
          std::forward<G> (gain));
        std::get<Idx> (*_delay).push_block (xspan {z.data(), in.size()});
      }
      else {
        fetch_block_optmod<Idx> (z, std::forward<Lfo> (lfo));
        std::get<Idx> (*_delay).push_block (in);
        xspan_memcpy (sout, z);
      }
    }
    else {
      // overlap, needs single-sample iteration. (Notice that it could be
      // smarter and run in the smallest possible block size, as of now not
      // worth the complexity, as most delays are bigger than the block size)
      auto lfo_gen = get_generic_generator (std::forward<Lfo> (lfo));
      auto g_gen   = get_gain_generator<Idx> (std::forward<G> (gain));
      constexpr std::optional<interpolation> interp = spec::get_interp (Idx);

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        value_type z;
        if constexpr (!interp) {
          z = read_next<Idx>();
        }
        else {
          // this is for completeness, modulations under the block size are
          // unlikely.
          z = in[i];
          if constexpr (*interp == interpolation::linear) {
            run_lerp<Idx> (xspan {&z, 1}, [i, &lfo_gen] (uint) {
              return lfo_gen (i);
            });
          }
          else {
            static_assert (*interp == interpolation::thiran);
            run_thiran<Idx> (xspan {&z, 1}, [i, &lfo_gen] (uint) {
              return lfo_gen (i);
            });
          }
        }
        if constexpr (spec::is_allpass (Idx) || spec::is_comb (Idx)) {
          if constexpr (spec::is_allpass (Idx)) {
            std::tie (out[i], z)
              = run_allpass<Idx, value_type> (in[i], z, g_gen (i));
          }
          else {
            std::tie (out[i], z)
              = run_comb<Idx, value_type> (in[i], z, g_gen (i));
          }
          std::get<Idx> (*_delay).push (z);
        }
        else {
          std::get<Idx> (*_delay).push (in[i]);
          out[i] = z;
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Tag, class... P>
  void run_multitap_delay (
    xspan<value_type const> in,
    Tag,
    value_type* out,
    P&&... outs)
  {
    constexpr bool out_overwrite = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_spls const> spls = spec::get_delays_spls (Idx);
    mp11::mp_for_each<mp11::mp_iota_c<spls.size()>> ([&] (auto i) {
      static_assert (spls[i].to_int() >= max_block_size);
    });
    // in and out can be aliased, but none of the outs shall be aliased with
    // in or within.
    block_array                                         tap0;
    std::array<value_type * artv_restrict, spls.size()> tap_ptrs = std::apply (
      [&] (auto&&... args) {
        return make_array<value_type * artv_restrict> (
          tap0.data(), static_cast<value_type*> (&args[0])...);
      },
      std::forward_as_tuple (outs...));

    assert (in.size() <= max_block_size);
    for (uint i = 0; i < tap_ptrs.size(); ++i) {
      uint spls_u = spls[i].to_int();
      if constexpr (out_overwrite) {
        auto dst = xspan {tap_ptrs[i], in.size()};
        assert (dst.data() != in.data()); // unwanted aliasing with input
        std::get<Idx> (*_delay).read_block (
          xspan {dst.data(), in.size()}, spls_u);
      }
      else {
        xspan dst {tap_ptrs[i], in.size()};
        assert (dst.data() != in.data()); // unwanted aliasing with input
        block_array tmp_mem;
        xspan       tmp {tmp_mem.data(), in.size()};
        std::get<Idx> (*_delay).read_block (
          xspan {tmp.data(), in.size()}, spls_u);
        span_add (dst, tmp);
      }
    }
    std::get<Idx> (*_delay).push_block (in);
    if constexpr (out_overwrite) {
      xspan_memdump (out, xspan {tap0.data(), in.size()});
    }
    else {
      span_add (out, xspan {tap0.data(), in.size()});
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Tag, class... Args>
  void run_parallel_delays (xspan<value_type const> in, Tag, Args&&... outs_arg)
  {
    constexpr bool out_overwrite     = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_t const> g = spec::get_gains (Idx);
    constexpr auto                 n_outs       = spec::get_n_outs (Idx);
    constexpr auto                 n_extra_buff = out_overwrite ? n_outs : 0;
    static_assert (n_outs <= g.size());
    static_assert (sizeof...(outs_arg) == n_outs);

    mp11::mp_repeat_c<std::tuple<block_array>, g.size() - n_extra_buff> tapmem;
    std::apply (
      [&] (auto&&... args) {
        if constexpr (out_overwrite) {
          run_multitap_delay<Idx> (
            in, overwrite, &outs_arg[0]..., args.data()...);
        }
        else {
          run_multitap_delay<Idx> (in, overwrite, args.data()...);
        }
      },
      tapmem);
    std::array<value_type*, g.size()> tap_ptrs = std::apply (
      [&] (auto&&... args) {
        if constexpr (out_overwrite) {
          return std::array {&outs_arg[0]..., args.data()...};
        }
        else {
          return std::array {args.data()...};
        }
      },
      tapmem);
    std::array<value_type*, n_outs> outs = std::apply (
      [&] (auto&&... args) { return std::array {&args[0]...}; },
      std::forward_as_tuple (outs_arg...));

    if constexpr (is_fixpt_v<value_type>) {
      // multiply and accumulate at high resolution
      std::array<std::array<fixpt_acum_t, max_block_size>, n_outs> acum {};
      for (uint tap = 0; tap < tap_ptrs.size(); ++tap) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          acum[tap % n_outs][i] += tap_ptrs[tap][i] * (value_type) g[tap];
        }
      }
      // quantize
      state::quantizer_arr_fixpt<n_outs>& st  = std::get<Idx> (_states);
      auto&                               err = st.err;

      for (uint out = 0; out < n_outs; ++out) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          auto [v, err_] = round (acum[out][i], err[out]);
          if constexpr (out_overwrite) {
            outs[out][i] = v;
          }
          else {
            outs[out][i] = (value_type) (v + outs[out][i]);
          }
          err[out] = err_;
        }
      }
    }
    else {
      for (uint tap = 0; tap < n_extra_buff; ++tap) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          // reminder: &out == &tap_ptrs[0]
          outs[tap][i] *= (value_type) g[tap];
        }
      }
      for (uint tap = n_extra_buff; tap < tap_ptrs.size(); ++tap) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          outs[tap % n_outs][i] += tap_ptrs[tap][i] * (value_type) g[tap];
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Range>
  void run_variable_delay (
    xspan<value_type>       out,
    xspan<value_type const> in,
    Range&&                 gen)
  {
    constexpr uint min   = spec::get_min_delay_spls (Idx);
    constexpr uint max   = spec::get_max_delay_spls (Idx);
    constexpr uint range = max - min;
    auto           rgen  = get_sample_generator (std::forward<Range> (gen));

    static_assert (min >= max_block_size);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      value_type v;
      uint       spls;
      auto       gv = rgen (i);
      if constexpr (spec::uses_absolute_delay_modulation (Idx)) {
        spls = (uint) gv;
      }
      else if constexpr (is_fixpt_v<value_type>) {
        spls = (uint) (min + range * gv.to_floatp());
      }
      else {
        spls = (uint) (min + range * gv);
      }
      assert (spls >= min);
      assert (spls <= max);
      v = get<Idx> (spls);
      std::get<Idx> (*_delay).push (in[i]);
      out[i] = v; // they could be aliased...
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Range>
  void run_variable_delay (xspan<value_type> io, Range&& gen)
  {
    run_variable_delay<Idx> (io, io, std::forward<Range> (gen));
  }
  //----------------------------------------------------------------------------
  template <std::size_t N>
  static constexpr int get_previous_ap_pos (
    uint                       pos,
    std::array<uint, N> const& idxs)
  {
    if (pos == 0) {
      return -1;
    }
    for (uint i = pos - 1; i < pos; --i) {
      if (spec::is_allpass (idxs[i])) {
        return i;
      }
    }
    return -1;
  }
  //----------------------------------------------------------------------------
  template <std::size_t N>
  static constexpr std::array<uint, N + 1> get_arg_offsets (
    std::array<uint, N> idxs)
  {
    uint                    arg = 0;
    std::array<uint, N + 1> ret {};
    for (uint i = 0; i < idxs.size(); ++i) {
      if (spec::is_allpass (idxs[i])) {
        // lfo + gain parameter.
        arg += 2;
        ret[i + 1] = arg;
        continue;
      }
      else if (spec::is_1tap_delay (idxs[i])) {
        // lfo for modulation
        arg += 1;
        ret[i + 1] = arg;
        continue;
      }
      else if (spec::is_lowpass (idxs[i]) || spec::is_highpass (idxs[i])) {
        // gain parameter.
        arg += 1;
        ret[i + 1] = arg;
        continue;
      }
      else if (spec::get_crossover_n_bands (idxs[i]) != 0) {
        // [freq, gain]... gain
        arg += (spec::get_crossover_n_bands (idxs[i]) * 2) + 1;
        ret[i + 1] = arg;
        continue;
      }
      else {
        assert (false); // Type not supported/added yet
      }
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  // run arbitrarily nested allpass.
  // For each stage all its arguments have to be given, those are
  // consumed positionally.
  //
  // Not provided parameters will be defaulted. Use "defaulted_tag" instances
  // to fill with defaults until the one you want to provide.
  //----------------------------------------------------------------------------
  template <uint... Idx, class... Ts>
  void run_nested_ap (
    value_type*             out,
    xspan<value_type const> in,
    Ts&&... argsp)
  {
    constexpr auto idxs        = make_array (Idx...);
    constexpr auto arg_offset  = get_arg_offsets (idxs);
    constexpr uint max_n_args  = arg_offset.back();
    constexpr int  last_ap_pos = get_previous_ap_pos (idxs.size(), idxs);

    static_assert (
      spec::is_allpass (idxs[0]), "1st stage has to be an allpass");
    static_assert (max_n_args >= sizeof...(Ts), "Excess arguments");
    using defaults = mp11::
      mp_repeat_c<std::tuple<defaulted_tag>, max_n_args - sizeof...(Ts)>;

    block_array fwd;
    auto args = std::tuple_cat (std::forward_as_tuple (argsp...), defaults {});
    xspan_memdump (fwd.data(), in);

    using idxlist = mp_list<k_uint<Idx>...>;

    mp_foreach_idx (idxlist {}, [&] (auto order, auto stage_idx) {
      bool constexpr is_serial = false; // not nested element // TBD

      if constexpr (spec::is_allpass (stage_idx) && !is_serial) {
        // having delays bigger than the block size allows to run each stage
        // once for all samples instead of each sample having to loop over
        // each stage.
        static_assert (
          spec::get_min_delay_spls (stage_idx) >= max_block_size,
          "Nested AP delays have to be GE than the block size");

        block_array bwd;
        fetch_block_optmod<stage_idx> (
          xspan {bwd.data(), in.size()}, std::get<arg_offset[order]> (args));

        auto&& gain_arg = std::get<arg_offset[order] + 1> (args);
        auto   gain_gen = get_gain_generator<stage_idx> (gain_arg);

        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          // See lattice form:
          // https://www.dsprelated.com/freebooks/pasp/Nested_Allpass_Filters.html
          std::tie (bwd[i], fwd[i])
            = run_allpass<stage_idx, value_type> (fwd[i], bwd[i], gain_gen (i));
        }
        if constexpr (order.value == 0) {
          // first allpass stage: backwards signal is the output
          xspan_memdump (out, xspan {bwd.data(), in.size()});
        }
        else {
          // try to process the backwards signal with nested stages in
          // backwards order until another allpass is found.
          constexpr int ap_pos = get_previous_ap_pos (order.value, idxs);
          static_assert (ap_pos >= 0);
          constexpr uint n_steps = (order - ap_pos) - 1;
          using seq = add_offset_t<ap_pos, std::make_index_sequence<n_steps>>;

          // process backward signal with nested non-allpass elements
          mp11::mp_for_each<mp11::mp_from_sequence<seq>> ([&] (auto order_rev) {
            constexpr uint rev_idx       = idxs[order_rev];
            constexpr bool is_serial_rev = false; // TBD
            if constexpr (is_serial_rev) {
              // skip... not belonging to the backwards path
            }
            else if constexpr (
              spec::is_highpass (rev_idx) || spec::is_lowpass (rev_idx)) {
              auto gain = std::get<arg_offset[order_rev]> (args);
              run (stage_list<rev_idx> {}, xspan {bwd.data(), in.size()}, gain);
            }
            else if constexpr (spec::is_1tap_delay (rev_idx)) {
              auto lfo_gen_rev = get_generic_generator (
                std::get<arg_offset[order_rev]> (args));
              run (
                stage_list<rev_idx> {},
                xspan {bwd.data(), in.size()},
                lfo_gen_rev);
            }
            else if constexpr (spec::get_crossover_n_bands (rev_idx) == 1) {
              static constexpr uint a1 = arg_offset[order_rev];
              run (
                stage_list<rev_idx> {},
                xspan {bwd.data(), in.size()},
                std::get<a1> (args),
                std::get<a1 + 1> (args),
                std::get<a1 + 2> (args));
            }
            else if constexpr (spec::get_crossover_n_bands (rev_idx) == 2) {
              static constexpr uint a1 = arg_offset[order_rev];
              run (
                stage_list<rev_idx> {},
                xspan {bwd.data(), in.size()},
                std::get<a1> (args),
                std::get<a1 + 1> (args),
                std::get<a1 + 2> (args),
                std::get<a1 + 3> (args),
                std::get<a1 + 4> (args));
            }
            else {
              static_assert (sizeof (order_rev), "Not implemented yet");
            }
          });
          // insert the backwards signal on the previous allpass...
          std::get<idxs[ap_pos]> (*_delay).push_block (
            xspan {bwd.data(), in.size()});
        }
      } // is_allpass
      else if constexpr ((order > last_ap_pos) || is_serial) {
        // elements after the last allpass can be processed on the forward
        // path so they are added to the last allpass on the final insertion
        // outside this foreach.
        //
        // If implementing serial elements inside the nestings they will
        // have to take this branch too, but instead of allpass they will
        // have to take a class like "serial<allpass> and be unwrapped"
        if constexpr (
          spec::is_highpass (stage_idx) || spec::is_lowpass (stage_idx)) {
          auto gain = std::get<arg_offset[order]> (args);
          run (stage_list<stage_idx> {}, xspan {fwd.data(), in.size()}, gain);
        }
        else if constexpr (spec::is_1tap_delay (stage_idx)) {
          auto lfo_gen_rev
            = get_generic_generator (std::get<arg_offset[order]> (args));
          run (
            stage_list<stage_idx> {},
            xspan {fwd.data(), in.size()},
            lfo_gen_rev);
        }
        else if constexpr (spec::get_crossover_n_bands (stage_idx) == 1) {
          static constexpr uint a1 = arg_offset[order];
          run (
            stage_list<stage_idx> {},
            xspan {fwd.data(), in.size()},
            std::get<a1> (args),
            std::get<a1 + 1> (args),
            std::get<a1 + 2> (args));
        }
        else if constexpr (spec::get_crossover_n_bands (stage_idx) == 2) {
          static constexpr uint a1 = arg_offset[order];
          run (
            stage_list<stage_idx> {},
            xspan {fwd.data(), in.size()},
            std::get<a1> (args),
            std::get<a1 + 1> (args),
            std::get<a1 + 2> (args),
            std::get<a1 + 3> (args),
            std::get<a1 + 4> (args));
        }
        else {
          static_assert (sizeof (stage_idx), "Not implemented yet");
        }
      }
      else {
        // skipping nested non-allpass elements, these are handdebugled inside
        // the. allpass conditional.
      }
    });
    // the last allpass gets the fwd signal enqueued
    static_assert (last_ap_pos > 0);
    std::get<idxs[last_ap_pos]> (*_delay).push_block (
      xspan {fwd.data(), in.size()});
  }
  //----------------------------------------------------------------------------
  // This might belong somewhere else if made more generic. Do when required.
  template <bool round, bool dither, class value_type>
  std::tuple<fixpt_t, s16> quantize (value_type spl, s16 err_prev)
  {
    static_assert (is_fixpt_v<value_type>);
    // Fraction-saving quantization
    // https://dsp.stackexchange.com/questions/66171/single-pole-iir-filter-fixed-point-design
    // https://dsp.stackexchange.com/questions/21792/best-implementation-of-a-real-time-fixed-point-iir-filter-with-constant-coeffic
    // noise shaping only
    constexpr uint n_truncated   = fixpt_acum_t::n_bits - fixpt_t::n_bits;
    constexpr uint mask          = lsb_mask<uint> (n_truncated);
    constexpr bool bidirectional = round || dither;

    if constexpr (dither) {
      spl = fixpt_clamp (
        spl,
        fixpt_t::min() + 3_r * fixpt_t::epsilon(),
        fixpt_t::max() - 3_r * fixpt_t::epsilon());
    }
    else {
      spl = fixpt_clamp (
        spl,
        fixpt_t::min() + fixpt_t::epsilon(),
        fixpt_t::max() - fixpt_t::epsilon());
    }

    auto         out = (fixpt_acum_t) (spl + fixpt_acum_t::from (err_prev));
    fixpt_acum_t d_out;
    if constexpr (dither) {
      // not throughly tested...
      auto noise = tpdf_dither<n_truncated, fixpt_acum_t::scalar_type> (_noise);
      d_out      = out + fixpt_acum_t::from (noise[0]);
    }
    else {
      d_out = out;
    }
    assert_range (d_out);
    fixpt_t q_out;
    if constexpr (round) {
      q_out = (fixpt_tr) d_out;
    }
    else {
      q_out = (fixpt_t) d_out; // truncation
    }
    if constexpr (bidirectional) {
      fixpt_acum_t errv = out - q_out;
      auto         err  = (s16) (errv.value() & mask);
      return std::make_tuple ((fixpt_t) q_out, err);
    }
    else {
      auto err = (s16) (out.value() & mask);
      return std::make_tuple (q_out, err);
    }
  }
  //----------------------------------------------------------------------------
  template <class value_type>
  std::tuple<fixpt_t, s16> round (value_type v, s16 err_prev)
  {
    return quantize<true, false> (v, err_prev);
  }
  //----------------------------------------------------------------------------
  template <class value_type>
  std::tuple<fixpt_t, s16> truncate (value_type v, s16 err_prev)
  {
    return quantize<false, false> (v, err_prev);
  }
  //----------------------------------------------------------------------------
  template <
    class value_type,
    std::enable_if_t<!is_fixpt_v<value_type>>* = nullptr>
  void assert_range (value_type v)
  {
#ifndef NDEBUG
    assert (abs (v) < 1.f);
#endif
  }
  //----------------------------------------------------------------------------
  template <
    class value_type,
    std::enable_if_t<is_fixpt_v<value_type>>* = nullptr>
  void assert_range (value_type v)
  {
#ifndef NDEBUG
    assert_range (v.to_floatp());
#endif
  }
  //----------------------------------------------------------------------------
  template <class Candidate>
  static constexpr bool is_generator()
  {
    using C = std::remove_reference_t<Candidate>;
    return std::is_convertible_v<C, std::function<value_type (uint)>>
      && !std::is_same_v<C, std::nullptr_t>;
  }
  //----------------------------------------------------------------------------
  template <class Candidate>
  static constexpr bool is_uint_generator()
  {
    using C = std::remove_reference_t<Candidate>;
    return std::is_convertible_v<C, std::function<uint (uint)>>
      && !std::is_same_v<C, std::nullptr_t>;
  }
  //----------------------------------------------------------------------------
  static constexpr void span_add (
    xspan<value_type>       dst,
    xspan<value_type const> lhs,
    xspan<value_type const> rhs)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < rhs.size(); ++i) {
      dst[i] = (value_type) (lhs[i] + rhs[i]);
    }
  }
  //----------------------------------------------------------------------------
  static constexpr void span_add (
    xspan<value_type>       lhs,
    xspan<value_type const> rhs)
  {
    span_add (lhs, lhs, rhs);
  }
  //----------------------------------------------------------------------------
  using block_array = std::array<value_type, max_block_size>;

  using indexes = mp11::mp_iota_c<spec::size()>;

  using states_typelist = std::conditional_t<
    is_fixpt_v<value_type>,
    mp11::mp_transform_q<state::index_to_state_qfn<fixpt_t, spec>, indexes>,
    mp11::mp_transform_q<state::index_to_state_qfn<float, spec>, indexes>>;
  using states_tuple = mp11::mp_rename<states_typelist, std::tuple>;

  using delays_typelist = mp11::
    mp_transform_q<delay::index_to_delay_buffer_qfn<Data_type, spec>, indexes>;
  using delays_tuple = mp11::mp_rename<delays_typelist, std::tuple>;

  delays_tuple*     _delay {nullptr};
  states_tuple      _states {};
  lowbias32_hash<1> _noise {};
};
}}} // namespace artv::detail::lofiverb
//------------------------------------------------------------------------------
