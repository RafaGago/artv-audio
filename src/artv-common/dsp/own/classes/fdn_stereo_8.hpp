#pragma once

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/diffusion_matrix.hpp"
#include "artv-common/dsp/own/classes/reverb_tools.hpp"

namespace artv {

//------------------------------------------------------------------------------
class fdn_stereo_8 {
public:
  //----------------------------------------------------------------------------
  static constexpr uint blocksize  = 16;
  static constexpr uint n_channels = 2;
  //----------------------------------------------------------------------------
  struct src_cfg {
    u16 kaiser_att_db;
    u16 taps_branch;
    u16 taps_branch_frac;
    u16 cutoff;
    u16 srate;
  };
  //----------------------------------------------------------------------------
  struct pre_dif_cfg {
    static constexpr uint              n_channels = 2;
    static constexpr uint              n_stages   = 4;
    float                              g_mod_depth;
    float                              g_mod_freq;
    float                              g_drift_fact;
    float                              g_max;
    array2d<u16, n_channels, n_stages> n_samples;
  };
  //----------------------------------------------------------------------------
  struct early_cfg {
    static constexpr uint n_channels = 4;
    static constexpr uint n_stages   = 4;

    u16   prime_idx;
    float rounding_factor;

    struct {
      float meters;
      float span;
      float g;
    } stage[n_stages];
  };
  //----------------------------------------------------------------------------
  struct late_cfg {
    static constexpr uint              n_channels = 16;
    static constexpr uint              n_stages   = 1;
    float                              rounding_factor;
    float                              span_factor;
    u16                                prime_idx;
    array2d<u16, n_channels, n_stages> n_samples;
  };
  //----------------------------------------------------------------------------
  struct int_dif_cfg {
    static constexpr uint              n_channels = 2;
    static constexpr uint              n_stages   = 4;
    float                              g_mod_depth;
    float                              g_mod_freq;
    float                              g_drift_fact;
    float                              g_base;
    u8                                 channel_l;
    u8                                 channel_r;
    array2d<u16, n_channels, n_stages> n_samples;
  };
  //----------------------------------------------------------------------------
  struct out_dif_cfg {
    static constexpr uint              n_channels = 2;
    static constexpr uint              n_stages   = 4;
    array2d<u16, n_channels, n_stages> n_samples;
  };
  //----------------------------------------------------------------------------
  struct cfg {
    src_cfg     src;
    pre_dif_cfg pre_dif;
    early_cfg   early;
    late_cfg    late;
    int_dif_cfg int_dif;
    out_dif_cfg out_dif;
  };
  //----------------------------------------------------------------------------
  static cfg get_default_cfg_preset()
  {
    // (1 + sqrt(5)) / 2;
    static constexpr double golden_ratio = 1.618033988749895;

    cfg r;
    r.src.kaiser_att_db    = 210;
    r.src.taps_branch      = 32;
    r.src.taps_branch_frac = 16;
    r.src.cutoff           = 9000;
    r.src.srate            = 27000;

    r.pre_dif.g_mod_depth  = 0.1f;
    r.pre_dif.g_mod_freq   = 1.6f;
    r.pre_dif.g_drift_fact = 0.8f;
    r.pre_dif.g_max        = 0.39f;

    r.pre_dif.n_samples = make_array (
      array_cast<u16> (make_array (43, 43)),
      array_cast<u16> (make_array (113, 113)),
      array_cast<u16> (make_array (317, 317)),
      array_cast<u16> (make_array (907, 906)));

    r.early.stage[0].meters = 2.f;
    r.early.stage[0].span   = 7.f;
    r.early.stage[0].g      = 0.9f;
    r.early.stage[1].meters = 2.4f;
    r.early.stage[1].span   = 3.f;
    r.early.stage[1].g      = 0.6f;
    r.early.stage[2].meters = 2.6f;
    r.early.stage[2].span   = 3.1f;
    r.early.stage[2].g      = 0.6f;
    r.early.stage[3].meters = 3.8f;
    r.early.stage[3].span   = 3.1f;
    r.early.stage[3].g      = 0.6f;
    r.early.prime_idx       = 1;
    r.early.rounding_factor = 1000;

    r.late.prime_idx       = 15;
    r.late.rounding_factor = 1;
    r.late.span_factor     = golden_ratio * 1.5f;
    r.late.n_samples       = make_array (array_cast<u16> (make_array (
      911,
      967,
      1181,
      1103,
      1289,
      1307,
      1669,
      1553,
      1753,
      1877,
      2131,
      2017,
      2647,
      2411,
      2957,
      2837)));

    r.int_dif.g_mod_depth  = 0.11f;
    r.int_dif.g_mod_freq   = 0.6f;
    r.int_dif.g_drift_fact = 0.8f;
    r.int_dif.n_samples    = make_array (
      array_cast<u16> (make_array (887, 887)),
      array_cast<u16> (make_array (478, 478)),
      array_cast<u16> (make_array (419, 421)),
      array_cast<u16> (make_array (907, 906)));

    r.out_dif.n_samples = make_array (
      array_cast<u16> (make_array (19, 23)),
      array_cast<u16> (make_array (53, 67)),
      array_cast<u16> (make_array (157, 191)),
      array_cast<u16> (make_array (443, 532)));

    return r;
  }
  //----------------------------------------------------------------------------
  void set_input_diffusor_gain (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _pre_dif_g = _cfg.pre_dif.g_max * factor;
  }
  //----------------------------------------------------------------------------
  void set_early_gain (float g) { _early_gain = g; }
  //----------------------------------------------------------------------------
  void reset (uint samplerate, cfg const& cfg)
  {
    _cfg = cfg;

    _resampler.reset (
      cfg.src.srate,
      samplerate,
      cfg.src.cutoff,
      cfg.src.cutoff,
      cfg.src.taps_branch,
      cfg.src.taps_branch_frac,
      cfg.src.kaiser_att_db,
      true,
      blocksize,
      6 * 1024);

    setup_early ((float) samplerate);
    memory_setup();
    // TODO: check no value under (blocksize) 16 after modulation
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    _resampler.process (outs, ins, samples, [=] (auto io) {
      process_block (io);
    });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void process_block (crange<std::array<float, 2>> io)
  {
    // hardcoding 2's to ease readability, as this will be always be stereo so
    // "static_asserting"
    static_assert (n_channels == 2);

    // the optimizer should be able to reduce the number of arrays used, as some
    // don't require persistency.
    array2d<float, 2, blocksize> tmp;
    array2d<float, 2, blocksize> pre_dif;
    array2d<float, 2, blocksize> early;
    static_assert (early_cfg::n_channels == 4); // hardcoding for readability
    array2d<float, 4, blocksize> early_mtx;
    array2d<float, 2, blocksize> late;

    while (io.size()) {
      auto block = io.cut_head (std::min<uint> (io.size(), blocksize));

      // pre diffusor ----------------------------------------------------------
      // AP gain LFO
      auto mod_g = make_crange (tmp);
      for (uint i = 0; i < block.size(); ++i) {
        vec<float, 2> mod = _pre_dif_lfo.tick_filt_sample_and_hold();
        vec<float, 2> g   = vec_set<2> (_pre_dif_g);
        g += mod * _cfg.pre_dif.g_mod_depth;
        mod_g[i] = vec_to_array (g);
      }
      // allpasses (todo loop ordering)
      for (uint i = 0; i < block.size(); ++i) {
        // hardcoding 2's to avoid unreadability, so "static_asserting"
        static_assert (n_channels == 2);
        pre_dif[i] = block[i];

        for (uint st = 0; st < _pre_dif.size(); ++st) {
          allpass_stage_tick (
            pre_dif[i], _pre_dif[st], mod_g[i], _cfg.pre_dif.n_samples[st]);
        }
      }
      // early -----------------------------------------------------------------
      // the stages are expanded blockwise manually as a perf optimization
      static_assert (early_cfg::n_stages == 4);

      // building the 4-wide matrices
      for (uint i = 0; i < block.size(); ++i) {
        early_mtx[i][0] = pre_dif[i][0];
        early_mtx[i][3] = pre_dif[i][1];
        early_mtx[i][1] = (early_mtx[i][0] + early_mtx[i][3]) * 0.5f; // mid
        early_mtx[i][2] = early_mtx[i][1];
      }

      for (uint stage = 0; stage < early_cfg::n_stages; ++stage) {
        // allpass
        auto g = array_broadcast<4> (_cfg.early.stage[stage].g);
        for (uint i = 0; i < block.size(); ++i) {
          allpass_stage_tick (
            early_mtx[i], _early[stage], g, _early_delay_spls[stage]);
        }
        // diffusion
        if (stage < 3) {
          for (uint i = 0; i < block.size(); ++i) {
            auto mtx = hadamard_matrix<4>::tick<float> (early_mtx[i]);
            if (stage == 0) {
              crange_memset<float> (early_mtx[i], 0);
            }
            for (uint j = 0; j < mtx.size(); ++j) {
              early_mtx[i][j] += mtx[j];
            }
          }
        }
        else {
          for (uint i = 0; i < block.size(); ++i) {
            auto mtx = householder_matrix<4>::tick<float> (early_mtx[i]);
            for (uint j = 0; j < mtx.size(); ++j) {
              early_mtx[i][j] += mtx[j];
            }
          }
        }

        uint n_rot = 0;
        switch (stage) {
        case 0:
          n_rot = 0;
          break;
        case 1:
          n_rot = 2;
          break;
        case 2:
          n_rot = 1;
          break;
        case 3:
          n_rot = 1;
          break;
        default:
          break;
        }
        // rotation
        for (uint i = 0; i < block.size(); ++i) {
          std::rotate (
            early_mtx[i].begin(),
            early_mtx[i].begin() + n_rot,
            early_mtx[i].end());
        }
      }
      for (uint i = 0; i < block.size(); ++i) {
        early[i][0] = early_mtx[i][0];
        early[i][1] = early_mtx[i][3];
      }

      // mixing ----------------------------------------------------------------
      for (uint i = 0; i < block.size(); ++i) {
        block[i][0] = early[i][0] * _early_gain;
        block[i][1] = early[i][1] * _early_gain;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void get_delay_length (
    crange<T> dst,
    T         spls_min,
    T         spls_max,
    uint      prime_idx,
    float     rounding_fact)
  {
    static constexpr uint n_max_stack = (6 * 1024) / sizeof (T);

    std::vector<T> dyn_mem;

    uint tbl_length = primes_table_size_guess (spls_min, spls_max);
    if (tbl_length > n_max_stack) {
      dyn_mem.resize (tbl_length);
      tbl_length = 0;
    }
    // Non portable: VLA.
    T stat_mem[tbl_length];

    auto work_mem = tbl_length ? make_crange (&stat_mem[0], tbl_length)
                               : make_crange (dyn_mem);

    delay_length::get (
      dst, spls_min, spls_max, prime_idx, rounding_fact, work_mem);
  }
  //----------------------------------------------------------------------------
  void setup_early (float srate)
  {
    // compute delay lengths
    for (uint i = 0; i < _cfg.early.n_stages; ++i) {
      uint spls_min = delay_length::meters_to_samples (
        _cfg.early.stage[i].meters, _cfg.src.srate, false);
      uint spls_max = delay_length::meters_to_samples (
        _cfg.early.stage[i].meters * _cfg.early.stage[i].span,
        _cfg.src.srate,
        true);
      get_delay_length<u16> (
        _early_delay_spls[i],
        spls_min,
        spls_max,
        _cfg.early.prime_idx,
        _cfg.early.rounding_factor);
    }
    // lfo
    _pre_dif_lfo.reset();
    auto f = vec_set<2> (_cfg.pre_dif.g_mod_freq);
    f[0] *= 1.01; // desync
    _pre_dif_lfo.set_freq (f, srate);
  };
  //----------------------------------------------------------------------------
  template <class T, size_t A, size_t B>
  static uint round_array_pow2_and_accumulate (array2d<T, A, B>& arr)
  {
    uint mem_total = 0;
    for (uint i = 0; i < B; ++i) {
      for (uint j = 0; j < A; ++j) {
        auto spls = arr[i][j];
        // check constraints for block processing, notice that this is not
        // accounting interpolation
        assert (spls > blocksize);
        auto new_spls = pow2_round_ceil (spls);
        assert (new_spls >= spls); // using uint16_t sometimes...
        arr[i][j] = new_spls;
        mem_total += new_spls;
      }
    }
    return mem_total;
  }
  //----------------------------------------------------------------------------
  void memory_setup()
  {
    auto pre_dif_n_spls = _cfg.pre_dif.n_samples;
    auto early_n_spls   = _early_delay_spls;

    // computing.
    uint mem_total = 0;
    mem_total += round_array_pow2_and_accumulate (pre_dif_n_spls);
    mem_total += round_array_pow2_and_accumulate (early_n_spls);

    // allocating
    _mem.clear();
    _mem.resize (mem_total);

    // assigning memory
    // pre diffusor
    auto mem = make_crange (_mem);
    for (uint i = 0; i < pre_dif_n_spls.size(); ++i) {
      for (uint j = 0; j < pre_dif_n_spls[0].size(); ++j) {
        _pre_dif[i][j].reset (mem.cut_head (pre_dif_n_spls[i][j]));
      }
    }
    // early
    for (uint i = 0; i < early_n_spls.size(); ++i) {
      for (uint j = 0; j < early_n_spls[0].size(); ++j) {
        _early[i][j].reset (mem.cut_head (early_n_spls[i][j]));
      }
    }
  }
  //----------------------------------------------------------------------------
  void allpass_stage_tick (
    crange<float>             io,
    crange<allpass<float_x1>> ap,
    crange<float>             g,
    crange<u16>               del_spls)
  {
    for (uint i = 0; i < io.size(); ++i) {
      io[i] = ap[i].tick (make_vec (io[i]), del_spls[i], make_vec (g[i]))[0];
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class Cfg>
  using reverb_array = array2d<T, Cfg::n_channels, Cfg::n_stages>;
  //----------------------------------------------------------------------------
  block_resampler<float, 2> _resampler;

  reverb_array<allpass<float_x1>, pre_dif_cfg> _pre_dif;
  lfo<2>                                       _pre_dif_lfo;
  float                                        _pre_dif_g;

  reverb_array<u16, early_cfg>               _early_delay_spls;
  reverb_array<allpass<float_x1>, early_cfg> _early;

  std::vector<float_x1> _mem;
  float                 _early_gain;

  cfg _cfg;
};
//------------------------------------------------------------------------------
} // namespace artv
