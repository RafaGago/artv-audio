#pragma once

// internals of lofiverb.hpp

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <vector>

#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/float.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

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
//------------------------------------------------------------------------------
struct delay_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g; // If 0 the delay has no allpass
};
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::delay_data make_dd (
  u16   spls,
  float g   = 0.,
  u16   mod = 0)
{
  assert (spls < fixpt_spls::max_int());
  assert (mod < fixpt_spls_mod::max_int());
  return delay_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g)};
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::delay_data make_ap (
  u16   spls,
  float g,
  u16   mod = 0)
{
  return make_dd (spls, g, mod);
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::delay_data make_delay (u16 spls, u16 mod = 0)
{
  return make_dd (spls, 0., mod);
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::delay_data make_damp (float g = 0)
{
  return make_dd (0, g);
}
//------------------------------------------------------------------------------
// A class to abstract 16-bit storage, queue access and common DSP operations
// when building reverbs based on allpass loops. Both on fixed and floating
// point.
template <class Spec, uint Max_block_size>
class engine {
public:
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = Max_block_size;
  //----------------------------------------------------------------------------
  static constexpr uint get_required_size()
  {
    uint ret = 0;
    for (uint i = 0; i < decltype (_stage) {}.size(); ++i) {
      ret += get_buffer_size (i);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<s16> mem)
  {
    for (uint i = 0; i < _stage.size(); ++i) {
      _stage[i].z = mem.cut_head (get_buffer_size (i)).data();
    }
  }
  //----------------------------------------------------------------------------
  // Pure delays with "size > blocksize" are placed before feedback loops to
  // enable block processing, so it is possible to fetch the future
  // feedbacks/outputs (minus the first) at once and then to push them all at
  // once. This is exactly what fetch/push accomplish.
  //
  // dst has to contain one element more at the head, where the previous output
  // will be placed. Once the feedback is applied to a current input, this head
  // feedback sample (last output) can be dropped.
  template <uint Idx, class T>
  void fetch_block (xspan<T> dst, uint negative_offset = 0)
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (dd.g.value() == 0, "Not possible on allpasses");
    static_assert (dd.mod.value() == 0, "Not possible on Modulated delays");
    assert (dst);
    decode_read (dst, get_read_buffers<Idx> (dst.size(), negative_offset));
  }
  // see comment on fetch_block_plus_one
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push (xspan<T const> src)
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (dd.g.value() == 0, "Not possible on allpasses");
    static_assert (dd.mod.value() == 0, "Not possible on Modulated delays");
    assert (src);
    encode_write (prepare_block_insertion<Idx> (src.size()), src);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_lp (xspan<T> io, U g)
  {
    assert (io);
    auto y1    = fetch_accumulator<Idx> (T {});
    using Acum = decltype (y1);
    for (uint i = 0; i < io.size(); ++i) {
      bool is_zero;
      if constexpr (is_fixpt_v<T>) {
        is_zero = (g.value() == 0);
      }
      else {
        is_zero = (g == 0.f);
      }
      auto gv = is_zero ? (Acum) get_gain<Idx, T>() : (Acum) g;
      auto y  = (Acum) io[i];
      y       = (1_r - gv) * io[i];
      y       = y + (y1 * gv);
      y1      = y;
      io[i]   = (T) v;
    }
    save_accumulator<Idx> (y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_lp (xspan<T> io)
  {
    run_lp<Idx> (io, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_hp (xspan<T> io, U g)
  {
    assert (io);
    auto y1    = fetch_accumulator<Idx> (T {});
    using Acum = decltype (y1);
    for (uint i = 0; i < io.size(); ++i) {
      bool is_zero;
      if constexpr (is_fixpt_v<T>) {
        is_zero = (g.value() == 0);
      }
      else {
        is_zero = (g == 0.f);
      }
      auto gv = is_zero ? (Acum) get_gain<Idx, T>() : (Acum) g;
      auto y  = ((Acum) io[i]) * gv;
      y += (1_r - gv) * y1;
      io[i] = (T) y;
      y1    = y;
    }
    save_accumulator<Idx> (y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_hp (xspan<T> io)
  {
    run_hp<Idx> (io, T {});
  }
  //----------------------------------------------------------------------------
  // run an allpass or pure delay.
  // Lfo and gain parameters can be either value generators (lambdas with T(int)
  // signature), spans or dummy types. When dummy types (no span and no lambda
  // types, e.g. nullptr) the parameters are defaulted (no delay modulation,
  // same allpass gain as in the specification)
  template <uint Idx, class T, class L, class G>
  void run (xspan<T> io, L lfo, G gain)
  {
    run_impl<Idx> (
      io,
      get_gain_generator<Idx, T> (std::forward<G> (gain)),
      get_lfo_generator<T> (std::forward<L> (lfo)));
  }
  //----------------------------------------------------------------------------
  // run an allpass or pure delay
  template <uint Idx, class T>
  void run (xspan<T> io)
  {
    run<Idx> (io, nullptr, nullptr);
  }
  //----------------------------------------------------------------------------
  // run 2 level nested allpass.
  // Lfo and gain parameters can be either value generators (lambdas with T(int)
  // signature), spans or dummy types. When dummy types (no span and no lambda
  // types, e.g. nullptr) the parameters are defaulted (no delay modulation,
  // same allpass gain as in the specification)
  template <
    uint Idx1,
    uint Idx2,
    class T,
    class L1,
    class G1,
    class L2,
    class G2>
  void run (xspan<T> io, L1 lfo1, G1 g1, L2 lfo2, G2 g2)
  {
    run_impl<Idx1, Idx2> (
      io,
      get_gain_generator<Idx1, T> (std::forward<G1> (g1)),
      get_lfo_generator<T> (std::forward<L1> (lfo1)),
      get_gain_generator<Idx2, T> (std::forward<G2> (g2)),
      get_lfo_generator<T> (std::forward<L2> (lfo2)));
  }
  //----------------------------------------------------------------------------
  // run a 2-level nested allpass
  template <uint Idx1, uint Idx2, class T>
  void run (xspan<T> io)
  {
    run<Idx1, Idx2> (io, nullptr, nullptr, nullptr, nullptr);
  }
  //----------------------------------------------------------------------------
  template <uint Start_idx, uint N, class T>
  void run_multichannel (std::array<T*, N> io, uint block_size)
  {
    mp11::mp_for_each<mp11::mp_iota_c<N>> ([&] (auto idx) {
      run<Start_idx + idx.value> (xspan {io[idx.value], block_size});
    });
  }
  //----------------------------------------------------------------------------
  template <class T, uint N>
  static void hadamard_no_norm (std::array<T*, N> io, uint block_size)
  {
    static_assert (is_pow2 (N));
    if constexpr (N > 1) {
      constexpr int half_n = N / 2;
      // all this superfluous array copying stuff should be removed by the
      // optimizer.
      hadamard (array_slice<0, half_n> (io), block_size);
      hadamard (array_slice<half_n, half_n> (io), block_size);

      for (uint n = 0; n < half_n; ++n) {
        uint n1 = n;
        uint n2 = n + half_n;
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < block_size; ++i) {
          T a       = io[n1][i];
          T b       = io[n2][i];
          io[n1][i] = (a + b);
          io[n2][i] = (a - b);
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T, uint N>
  static void hadamard (std::array<T*, N> io, uint block_size)
  {
    hadamard_no_norm (io, block_size);
    // normalize by 1/sqrtn, not using GCEM because this benefits from fixed
    // point constant literals.
    if constexpr (N == 2) {
      for (auto iov : io) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < block_size; ++i) {
          io[i] = (T) (io[i] * get_hadamard_ratio<N>());
        }
      }
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T, class Lfo>
  static constexpr auto get_lfo_generator (Lfo&& v)
  {
    if constexpr (
      std::is_convertible_v<Lfo, std::function<T (uint)>>
      && !std::is_same_v<Lfo, std::nullptr_t>) {
      return std::forward<Lfo> (v);
    }
    else if constexpr (is_xspan<Lfo>) {
      return [v] (uint i) { return v[i]; };
    }
    else {
      return nullptr;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class Gain>
  static constexpr auto get_gain_generator (Gain&& v)
  {
    if constexpr (
      std::is_convertible_v<Gain, std::function<T (uint)>>
      && !std::is_same_v<Gain, std::nullptr_t>) {
      return std::forward<Gain> (v);
    }
    else if constexpr (is_xspan<Gain>) {
      return [v] (uint i) { return v[i]; };
    }
    else {
      return [] (uint) { return get_gain<Idx, T>(); };
    }
  }
  //----------------------------------------------------------------------------
  template <uint N>
  static constexpr auto get_hadamard_ratio()
  {
    if constexpr (N == 4) {
      return 0.5_r;
    }
    else if constexpr (N == 8) {
      return 0.35355339059327373_r;
    }
    else if constexpr (N == 16) {
      return 0.25_r;
    }
    else if constexpr (N == 32) {
      return 0.17677669529663687_r;
    }
    else if constexpr (N == 64) {
      return 0.125_r;
    }
    else {
      static_assert (N == N, "Unimplemented");
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  static constexpr auto get_gain()
  {
    constexpr delay_data dd = Spec::values[Idx];
    if constexpr (is_fixpt_v<T>) {
      return dd.g;
    }
    else {
      return (float) dd.g;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void decode_read (T& dst, s16 src)
  {
    if constexpr (std::is_same_v<T, fixpt_t>) {
      dst = fixpt_sto::from (src);
    }
    else {
      static_assert (std::is_same_v<T, float>);
      dst = float16::decode (src);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void decode_read (xspan<T> dst, std::array<xspan<s16>, 2> src)
  {
    assert (dst.size() == (src[0].size() + src[1].size()));
    for (uint i = 0; i < src[0].size(); ++i) {
      decode_read (dst[i], src[0][i]);
    }
    dst.cut_head (src[0].size());
    for (uint i = 0; i < src[1].size(); ++i) {
      decode_read (dst[i], src[1][i]);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void encode_write (s16& dst, T src)
  {
    if constexpr (std::is_same_v<T, fixpt_t>) {
      dst = ((fixpt_sto) src).value();
    }
    else {
      static_assert (std::is_same_v<T, float>);
      dst = float16::encode (src);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void encode_write (std::array<xspan<s16>, 2> dst, xspan<T const> src)
  {
    assert (src.size() == (dst[0].size() + dst[1].size()));
    for (uint i = 0; i < dst[0].size(); ++i) {
      encode_write (dst[0][i], src[i]);
    }
    src.cut_head (dst[0].size());
    for (uint i = 0; i < dst[1].size(); ++i) {
      encode_write (dst[1][i], src[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void advance_pos()
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    ++_stage[Idx].pos;
    if (_stage[Idx].pos == sz) {
      _stage[Idx].pos = 0;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  T read_next()
  {
    constexpr delay_data dd = Spec::values[Idx];
    return get<Idx, T> (dd.spls.to_int());
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push_one (T v)
  {
    encode_write (_stage[Idx].z[_stage[Idx].pos], v);
    advance_pos<Idx>();
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  fixpt_acum_t fetch_accumulator (fixpt_t)
  {
    static_assert (sizeof (fixpt_acum_t) == 2 * sizeof (fixpt_t));
    fixpt_acum_t   ret;
    constexpr auto tail = get_delay_size (Idx);
    memcpy (&ret, &_stage[Idx].z[tail], sizeof ret);
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  float fetch_accumulator (float)
  {
    static_assert (sizeof (float) == 2 * sizeof (fixpt_t));
    float          ret;
    constexpr auto tail = get_delay_size (Idx);
    memcpy (&ret, &_stage[Idx].z[tail], sizeof ret);
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void save_accumulator (fixpt_acum_t v)
  {
    static_assert (sizeof (v) == 2 * sizeof (fixpt_t));
    constexpr auto tail = get_delay_size (Idx);
    memcpy (&_stage[Idx].z[tail], &v, sizeof v);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void save_accumulator (float v)
  {
    static_assert (sizeof (v) == 2 * sizeof (fixpt_t));
    constexpr auto tail = get_delay_size (Idx);
    memcpy (&_stage[Idx].z[tail], &v, sizeof v);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_thiran (xspan<fixpt_t> dst, LF&& lfo_gen)
  {
    constexpr delay_data dd     = Spec::values[Idx];
    constexpr auto       sz     = get_delay_size (Idx);
    constexpr auto       y1_pos = sz;

    fixpt_t y1;
    decode_read (y1, _stage[Idx].z[y1_pos]);
    for (uint i = 0; i < dst.size(); ++i) {
      auto fixpt_spls
        = (dd.spls.add_sign() + (lfo_gen (i) * dd.mod.add_sign()));
      auto n_spls = fixpt_spls.to_int();
      n_spls -= i;
      fixpt_t n_spls_frac = (fixpt_t) fixpt_spls.fractional();

      auto    d = (n_spls_frac + 0.418_r) & fixpt_resize_token<0, -1> {};
      fixpt_t a = (fixpt_t) ((1_r - d) / (1_r + d)); // 0.4104 to -1

      auto z0 = get<Idx, fixpt_t> (n_spls - 1);
      auto z1 = get<Idx, fixpt_t> (n_spls);
      y1      = (fixpt_t) (z0 * a + z1 - a * y1);
      dst[i]  = y1;
    }
    encode_write (_stage[Idx].z[y1_pos], y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_thiran (xspan<float> dst, LF&& lfo_gen)
  {
    constexpr delay_data dd     = Spec::values[Idx];
    constexpr auto       sz     = get_delay_size (Idx);
    constexpr auto       y1_pos = sz;

    float y1;
    decode_read (y1, _stage[Idx].z[y1_pos]);
    for (uint i = 0; i < dst.size(); ++i) {
      float fixpt_spls
        = dd.spls.to_floatp() + (lfo_gen (i) * dd.mod.to_floatp());
      uint  n_spls      = (uint) fixpt_spls;
      float n_spls_frac = fixpt_spls - n_spls;
      n_spls -= i;

      float d = n_spls_frac + 0.418_r;
      float a = (1_r - d) / (1_r + d); // 0.4104 to -1

      auto z0 = get<Idx, float> (n_spls - 1);
      auto z1 = get<Idx, float> (n_spls);
      y1      = z0 * a + z1 - a * y1;
      dst[i]  = y1;
    }
    encode_write (_stage[Idx].z[y1_pos], y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  T get (uint delay_spls)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    assert (delay_spls < sz);

    uint z = _stage[Idx].pos - delay_spls;
    z += (z >= sz) ? sz : 0;
    T ret;
    decode_read (ret, _stage[Idx].z[z]);
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  auto run_allpass (T in, T yn, U g)
  {
    auto u = in + yn * g;
    if constexpr (std::is_floating_point_v<T>) {
      if (abs (u) >= 1.) {
        assert (false);
      }
    }
    auto x = yn - u * g;
    if constexpr (std::is_floating_point_v<T>) {
      if (abs (x) >= 1.) {
        assert (false);
      }
    }
    return std::make_tuple (static_cast<T> (x), static_cast<T> (u));
  }
  //----------------------------------------------------------------------------
  // Template function to run:
  // - delays
  // - modulated delays
  // - non-nested allpasses
  // - modulated non-nested allpasses.
  //
  template <uint Idx, class T, class GF, class LF>
  void run_impl (xspan<T> io, GF&& g_gen, LF&& lfo_gen)
  {
    constexpr delay_data dd             = Spec::values[Idx];
    constexpr auto       sz             = get_delay_size (Idx);
    constexpr auto       minsz          = sz - dd.mod.to_int();
    constexpr auto       has_allpass    = dd.g.value() != 0;
    constexpr auto       has_modulation = dd.mod.value() != 0;

    auto szw = get_delay_size (Idx);

    if constexpr (minsz >= max_block_size) {
      // no overlap, can run block-wise
      assert (io.size() <= minsz);

      std::array<T, max_block_size> iocp_mem;
      xspan                         iocp {iocp_mem.data(), io.size()};
      xspan_memcpy (iocp, io);

      if constexpr (has_modulation) {
        run_thiran<Idx> (io, std::forward<LF> (lfo_gen));
      }
      else {
        decode_read (io, get_read_buffers<Idx> (io.size(), 0));
      }
      if constexpr (has_allpass) {
        for (uint i = 0; i < io.size(); ++i) {
          auto [out, push] = run_allpass<Idx, T> (iocp[i], io[i], g_gen (i));
          io[i]            = out;
          iocp[i]          = push;
        }
      }
      encode_write (prepare_block_insertion<Idx> (io.size()), iocp.to_const());
    }
    else {
      // overlap, needs single-sample iteration
      for (uint i = 0; i < io.size(); ++i) {
        T delayed;
        if constexpr (has_modulation) {
          // this is for completeness, modulations under the block size are
          // unlikely.
          delayed = io[i];
          run_thiran<Idx> (xspan {&delayed, 1}, [i, &lfo_gen] (uint) {
            return lfo_gen (i);
          });
        }
        else {
          delayed = read_next<Idx, T>();
        }
        if constexpr (has_allpass) {
          auto [out, push] = run_allpass<Idx, T> (io[i], delayed, g_gen (i));
          push_one<Idx> (push);
          io[i] = out;
        }
        else {
          push_one<Idx> (io[i]);
          io[i] = delayed;
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  // Helper template function to run (arbitrarily) nested allpasses
  template <
    uint I,
    uint N,
    class Idxs,
    class T,
    class GainGen,
    class LfoGen,
    class... GainAndLfoGens>
  void run_impl_get_yn (
    std::array<xspan<T>, N>& yn,
    GainGen&&                _,
    LfoGen&&                 lfo_gen,
    GainAndLfoGens&&... fwd)
  {
    constexpr uint       idx            = mp11::mp_at_c<Idxs, I>::value;
    constexpr delay_data dd             = Spec::values[idx];
    constexpr uint       sz             = get_delay_size (idx);
    constexpr uint       minsz          = sz - dd.mod.to_int();
    constexpr bool       has_modulation = dd.mod.value() != 0;

    // only block processing for now, non-block processing could be added.
    static_assert (minsz >= max_block_size);
    assert (yn.size() <= minsz);

    if constexpr (has_modulation) {
      run_thiran<idx> (yn[I], std::forward<LfoGen> (lfo_gen));
    }
    else {
      decode_read (yn[I], get_read_buffers<idx> (yn[I].size(), 0));
    }
    // recurse....
    if constexpr ((I + 1) < N) {
      run_impl_get_yn<I + 1, N, Idxs> (
        yn, std::forward<GainAndLfoGens> (fwd)...);
    }
  }
  //----------------------------------------------------------------------------
  // Helper template function to run nested allpasses
  template <
    uint I,
    uint N,
    class T,
    class GainGen,
    class LfoGen,
    class... GainAndLfoGens>
  void run_impl_get_gains (
    std::array<T, N>& gains,
    uint              n_spl,
    GainGen&&         gain_gen,
    LfoGen&&          _,
    GainAndLfoGens&&... fwd)
  {
    gains[I] = gain_gen (n_spl);
    // recurse....
    if constexpr ((I + 1) < N) {
      run_impl_get_gains<I + 1, N> (
        gains, n_spl, std::forward<GainAndLfoGens> (fwd)...);
    }
  }
  //----------------------------------------------------------------------------
  // Helper template function to run nested allpasses
  template <class T, uint N, class... GainAndLfoGens>
  std::array<T, N> run_impl_get_gains (uint n_spl, GainAndLfoGens&&... glfo)
  {
    std::array<T, N> ret;
    run_impl_get_gains<0, N> (
      ret, n_spl, std::forward<GainAndLfoGens> (glfo)...);
    return ret;
  }
  //----------------------------------------------------------------------------
  // Template function to run arbitrarily nested allpass/plain delay
  // combinations with or gain and delay time modulation (the gain and time
  // modulations can be unity but have to be present). It only work for delay
  // sizes bigger than the block size for implementation "simplicity".
  template <uint... Idx, class T, class... GainAndLfoGens>
  void run_impl (xspan<T> io, GainAndLfoGens&&... gnlfo)
  {
    using Idxs       = mp_list<std::integral_constant<uint, Idx>...>;
    constexpr uint n = mp11::mp_size<Idxs>::value;
    static_assert ((n * 2) == mp11::mp_size<mp_list<GainAndLfoGens...>>::value);

    array2d<T, max_block_size, n> yn_mem;
    std::array<xspan<T>, n>       yn;
    for (uint i = 0; i < n; ++i) {
      yn[i] = xspan {yn_mem[i].data(), io.size()};
    }
    run_impl_get_yn<0, n, Idxs> (yn, std::forward<GainAndLfoGens> (gnlfo)...);

    for (uint i = 0; i < io.size(); ++i) {
      auto g
        = run_impl_get_gains<T, n> (i, std::forward<GainAndLfoGens> (gnlfo)...);

      std::array<T, n> u;

      mp_foreach_idx (Idxs {}, [&] (auto j, auto topo_stage) {
        if (j == 0) {
          // The first stage will be always an allpass
          u[0]  = (T) (io[i] + yn[0][i] * g[0]);
          io[i] = (T) (yn[0][i] - u[0] * g[0]); // output
        }
        else {
          // subsequent stages of the lattice can be allpasses or plain delays
          if constexpr (Spec::values[topo_stage.value].mod.to_int() != 0) {
            u[0] = (T) (u[0] + yn[j][i] * g[j]);
            u[j] = (T) (yn[j][i] - u[0] * g[j]);
          }
          else {
            // The optimizer would see the 0 mul anyways, just for highlighting
            // that plain delays can be added to the allpass chain
            u[j] = yn[j][i];
          }
        }
      });
      for (uint j = 0; j < (n - 1); ++j) {
        yn[j][i] = u[j + 1];
      }
      yn[n - 1][i] = u[0];
    }
    mp_foreach_idx (Idxs {}, [&] (auto i, auto topo_stage) {
      encode_write (
        prepare_block_insertion<topo_stage.value> (io.size()),
        yn[i].to_const());
    });
  }
#if 0
  //----------------------------------------------------------------------------
  // Template function to run double nested allpasses, kept for reference
  // understanding what is going on on the implementation for arbitrary nested
  // allpasses above, as it is relatively "clever" code.
  template <
    uint Idx1,
    uint Idx2,
    class T,
    class GF1,
    class LF1,
    class GF2,
    class LF2>
  void run_impl (
    xspan<T> io,
    GF1&&    g_gen1,
    LF1&&    lfo_gen1,
    GF2&&    g_gen2,
    LF2&&    lfo_gen2)
  {
    constexpr delay_data dd[2] = {Spec::values[Idx1], Spec::values[Idx2]};
    constexpr uint       sz[2] = {get_delay_size (Idx1), get_delay_size (Idx2)};
    constexpr uint       minsz[2]
      = {sz[0] - dd[0].mod.to_int(), sz[1] - dd[1].mod.to_int()};
    constexpr bool has_modulation[2]
      = {dd[0].mod.value() != 0, dd[1].mod.value() != 0};

    // only block processing for now, non-block processing can be easily added.
    static_assert (minsz[0] >= max_block_size && minsz[1] >= max_block_size);

    assert ((io.size() <= minsz[0]) && (io.size() <= minsz[1]));

    array2d<T, max_block_size, 2> yn_mem;
    auto                          yn = make_array (
      xspan {yn_mem[0].data(), io.size()}, xspan {yn_mem[1].data(), io.size()});

    if constexpr (has_modulation[0]) {
      run_thiran<Idx1> (yn[0], std::forward<LF1> (lfo_gen1));
    }
    else {
      decode_read (yn[0], get_read_buffers<Idx1> (yn[0].size(), 0));
    }
    if constexpr (has_modulation[1]) {
      run_thiran<Idx2> (yn[1], std::forward<LF2> (lfo_gen2));
    }
    else {
      decode_read (yn[1], get_read_buffers<Idx2> (yn[1].size(), 0));
    }
    for (uint i = 0; i < io.size(); ++i) {
      auto g = make_array (g_gen1 (i), g_gen2 (i));

      auto u   = (T) (io[i] + yn[0][i] * g[0]);
      io[i]    = (T) (yn[0][i] - u * g[1]); // output
      u        = (T) (u + yn[1][i] * g[1]);
      auto u1  = (T) (yn[1][i] - u * g[1]);
      yn[0][i] = u1;
      yn[1][i] = u;
    }
    encode_write (prepare_block_insertion<Idx1> (io.size()), yn[0].to_const());
    encode_write (prepare_block_insertion<Idx2> (io.size()), yn[1].to_const());
  }
  //----------------------------------------------------------------------------
#endif
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_size (uint i)
  {
    uint ret = Spec::values[i].spls.to_int() + 1; // spls + 1 = size
    ret += Spec::values[i].mod.to_int(); // add mod spls to the size
    return ret;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_buffer_size (uint i)
  {
    uint ret = get_delay_size (i);
    ret += ((uint) (Spec::values[i].mod.value() != 0)) * 2; // thiran state
    ret += ((uint) (Spec::values[i].spls.value() == 0)) * 2; // damp filt state
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  std::array<xspan<s16>, 2> get_read_buffers (uint blocksize, uint neg_offset)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);
    assert (blocksize);
    assert (sz >= (blocksize + neg_offset));

    uint block1 = _stage[Idx].pos - dd.spls.to_int() - neg_offset;
    block1 += (block1 >= sz) ? sz : 0;
    uint end = _stage[Idx].pos - dd.spls.to_int() + blocksize - neg_offset - 1;
    end += (end >= sz) ? sz : 0;

    uint block1sz, block2sz;
    if (block1 < end) {
      // contiguous
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = sz - block1;
      block2sz = blocksize - block1sz;
    }
    return {
      xspan {&_stage[Idx].z[block1], block1sz},
      xspan {&_stage[Idx].z[0], block2sz}};
  }
  //----------------------------------------------------------------------------
  // moves the write pointer and returns the locations to write
  template <uint Idx>
  std::array<xspan<s16>, 2> prepare_block_insertion (uint blocksize)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);
    assert (blocksize);
    assert (sz >= blocksize);

    uint block1 = _stage[Idx].pos;
    uint end    = block1 + blocksize - 1;
    end -= (end >= sz) ? sz : 0;
    _stage[Idx].pos = end;
    advance_pos<Idx>();

    uint block1sz, block2sz;
    if (block1 < end) {
      // contiguous
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = sz - block1;
      block2sz = blocksize - block1sz;
    }
    return {
      xspan {&_stage[Idx].z[block1], block1sz},
      xspan {&_stage[Idx].z[0], block2sz}};
  }
  //----------------------------------------------------------------------------
  using float16 = f16pack<5, -1, f16pack_dftz>;
  struct stage {
    s16* z {};
    uint pos {};
  };
  std::array<stage, Spec::values.size()> _stage;
};
}}} // namespace artv::detail::lofiverb
//------------------------------------------------------------------------------
