#pragma once

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/diffusion_matrix.hpp"
#include "artv-common/dsp/own/classes/reverb_tools.hpp"
#include "artv-common/dsp/own/parts/interpolation/stateless.hpp"

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
    static constexpr uint       n_channels = 16;
    static constexpr uint       n_stages   = 1;
    float                       rounding_factor;
    float                       span_factor;
    float                       max_chorus_width;
    float                       max_chorus_freq;
    u16                         max_chorus_depth_spls;
    u16                         prime_idx;
    std::array<u16, n_channels> n_samples;
  };
  //----------------------------------------------------------------------------
  struct internal_dif_cfg {
    static constexpr uint              n_channels = 2;
    static constexpr uint              n_stages   = 4;
    float                              g_mod_depth;
    float                              g_mod_freq;
    float                              g_base;
    u8                                 channel_l;
    u8                                 channel_r;
    array2d<u16, n_channels, n_stages> n_samples;
  };
  //----------------------------------------------------------------------------
  struct out_dif_cfg {
    static constexpr uint              n_channels = 2;
    static constexpr uint              n_stages   = 4;
    float                              g_max;
    array2d<u16, n_channels, n_stages> n_samples;
  };
  //----------------------------------------------------------------------------
  struct filter_cfg {
    static constexpr uint       n_channels = 16;
    float                       max_att_db;
    float                       freq_factor;
    std::array<u16, n_channels> freqs;
  };
  //----------------------------------------------------------------------------
  struct cfg {
    src_cfg          src;
    pre_dif_cfg      pre_dif;
    early_cfg        early;
    late_cfg         late;
    internal_dif_cfg int_dif;
    out_dif_cfg      out_dif;
    filter_cfg       filter;
  };
  //----------------------------------------------------------------------------
  // a helper to aid in writing "cfg::late::n_samples" in ascending order.
  template <class T>
  static void from_ascending_pairs_to_internal_chnl_order (
    std::array<T, late_cfg::n_channels>& arr)
  {
    auto tmp = arr;

    arr[0]  = tmp[15];
    arr[1]  = tmp[13];
    arr[2]  = tmp[11];
    arr[3]  = tmp[9];
    arr[4]  = tmp[7];
    arr[5]  = tmp[5];
    arr[6]  = tmp[3];
    arr[7]  = tmp[1];
    arr[8]  = tmp[0];
    arr[9]  = tmp[2];
    arr[10] = tmp[4];
    arr[11] = tmp[6];
    arr[12] = tmp[8];
    arr[13] = tmp[10];
    arr[14] = tmp[12];
    arr[15] = tmp[14];
  }
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

    r.pre_dif.g_mod_depth = 0.1f;
    r.pre_dif.g_mod_freq  = 1.6f;
    r.pre_dif.g_max       = 0.39f;

    r.pre_dif.n_samples = make_array (
      array_cast<u16> (make_array (43, 43)),
      array_cast<u16> (make_array (113, 113)),
      array_cast<u16> (make_array (317, 317)),
      array_cast<u16> (make_array (907, 906)));

    r.early.stage[0].meters = 3.1f;
    r.early.stage[0].span   = 3.2f;
    r.early.stage[0].g      = 0.5f;
    r.early.stage[1].meters = 4.1f;
    r.early.stage[1].span   = 5.f;
    r.early.stage[1].g      = 0.3f;
    r.early.stage[2].meters = 2.1f;
    r.early.stage[2].span   = 3.1f;
    r.early.stage[2].g      = 0.35f;
    r.early.stage[3].meters = 1.9f;
    r.early.stage[3].span   = 3.1f;
    r.early.stage[3].g      = 0.8f;
    r.early.prime_idx       = 1;
    r.early.rounding_factor = 1000;

    r.late.prime_idx       = 15;
    r.late.rounding_factor = 1;
    r.late.span_factor     = golden_ratio * 1.5f;

    r.late.max_chorus_depth_spls = 150;
    r.late.max_chorus_freq       = 7;
    r.late.max_chorus_width      = 0.5f;

    r.late.n_samples = array_cast<u16> (make_array (
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
      2837));

    from_ascending_pairs_to_internal_chnl_order (r.late.n_samples);

    r.int_dif.g_mod_depth = 0.11f;
    r.int_dif.g_mod_freq  = 0.6f;
    r.int_dif.g_base      = 0.45f;
    r.int_dif.channel_l   = 2;
    r.int_dif.channel_r   = 13;

    r.int_dif.n_samples = make_array (
      array_cast<u16> (make_array (887, 887)),
      array_cast<u16> (make_array (478, 478)),
      array_cast<u16> (make_array (419, 421)),
      array_cast<u16> (make_array (907, 906)));

    r.out_dif.g_max     = 0.77f;
    r.out_dif.n_samples = make_array (
      array_cast<u16> (make_array (19, 23)),
      array_cast<u16> (make_array (53, 67)),
      array_cast<u16> (make_array (157, 191)),
      array_cast<u16> (make_array (443, 532)));

    r.filter.max_att_db  = -9.f;
    r.filter.freq_factor = -1.6f;

    r.filter.freqs = array_cast<u16> (make_array (
      500,
      530,
      640,
      600,
      860,
      800,
      920,
      900,
      1000,
      1100,
      1200,
      1000,
      2200,
      2000,
      2500,
      2700));

    from_ascending_pairs_to_internal_chnl_order (r.filter.freqs);

    return r;
  }
  //----------------------------------------------------------------------------
  void set_input_diffusor_gain (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _pre_dif_g = _cfg.pre_dif.g_max * factor;
  }
  //----------------------------------------------------------------------------
  void set_output_diffusor_gain (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _out_dif_g = _cfg.out_dif.g_max * factor;
  }
  //----------------------------------------------------------------------------
  void set_early_gain (float g) { _early_gain = g; }
  //----------------------------------------------------------------------------
  void set_late_gain (float g) { _late_gain = g; }
  //----------------------------------------------------------------------------
  void set_mod_freq (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _mod_freq_hz = factor * factor * _cfg.late.max_chorus_freq;
    reset_late_lfo();
  }
  //----------------------------------------------------------------------------
  void set_mod_depth (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _mod_depth_spls = factor * factor * (float) _cfg.late.max_chorus_depth_spls;
  }
  //----------------------------------------------------------------------------
  void set_mod_stereo (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    _mod_stereo = factor;
    reset_late_lfo();
  }
  //----------------------------------------------------------------------------
  void set_early_to_late (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _er_2_late = factor;
  }
  //----------------------------------------------------------------------------
  void set_in_to_late (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    _in_2_late = factor;
  }
  //----------------------------------------------------------------------------
  void set_l_matrix_angle (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _late_l_angle = set_rotation_angle (factor);
  }
  //----------------------------------------------------------------------------
  void set_r_matrix_angle (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _late_r_angle = set_rotation_angle (factor);
  }
  //----------------------------------------------------------------------------
  void set_lr_matrix_angle (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _late_lr_angle = set_rotation_angle (factor);
  }
  //----------------------------------------------------------------------------
  void set_time_msec (float msec)
  {
    _seconds = msec * 0.001f;
    reset_times();
  }
  //----------------------------------------------------------------------------
  void set_damp_freq (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    _damp_freq = -factor;
    reset_damping_filters();
  }
  //----------------------------------------------------------------------------
  void set_damp_factor (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _damp_factor = factor;
    reset_damping_filters();
  }
  //----------------------------------------------------------------------------
  void set_hp_freq (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    float hp_freq = 2.f + (factor * 58.f); // TODO: make variables
    _filters.reset_coeffs<dc_idx> (
      vec_set<16> (hp_freq), (float) _cfg.src.srate);
  }
  //----------------------------------------------------------------------------
  void set_lf_time_factor (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _lf_rt60_factor = 0.3333f + factor * 1.6666f; // TODO: make variables
    reset_times();
  }
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

    setup_in_diffusor();
    setup_early();
    setup_late();
    setup_internal_diffusor();
    setup_memory();

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

    static_assert (late_cfg::n_channels == 16); // hardcoding for readability
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
      for (uint i = 0; i < block.size(); ++i) {
        // hardcoding 2's to avoid unreadability, so "static_asserting"
        static_assert (n_channels == 2);
        pre_dif[i] = block[i];
      }

      for (uint st = 0; st < _pre_dif.size(); ++st) {
        for (uint i = 0; i < block.size(); ++i) {
          allpass_stage_tick (
            pre_dif[i], _pre_dif[st], mod_g[i], _cfg.pre_dif.n_samples[st]);
        }
      }

      // early -----------------------------------------------------------------
      // the stages are expanded blockwise manually as a perf optimization
      static_assert (early_cfg::n_stages == 4);

      memset (&early, 0, sizeof early);

      // building the 4-wide matrices
      for (uint i = 0; i < block.size(); ++i) {
        early_mtx[i][0] = pre_dif[i][0];
        early_mtx[i][1] = pre_dif[i][1];
        early_mtx[i][2] = (early_mtx[i][0] + early_mtx[i][1]) * 0.5f; // mid
        early_mtx[i][3] = early_mtx[i][2];
      }

      for (uint stage = 0; stage < 1 /*early_cfg::n_stages*/; ++stage) {
        // allpass
        auto g = array_broadcast<4> (_cfg.early.stage[stage].g);
        for (uint i = 0; i < block.size(); ++i) {
          allpass_stage_tick (
            early_mtx[i], _early[stage], g, _early_delay_spls[stage]);
        }
        // diffusion
        if (stage < 3) {
          for (uint i = 0; i < block.size(); ++i) {
            early_mtx[i] = hadamard_matrix<4>::tick<float> (early_mtx[i]);
          }
        }
        else {
          for (uint i = 0; i < block.size(); ++i) {
            early_mtx[i] = householder_matrix<4>::tick<float> (early_mtx[i]);
          }
        }
        // rotation
        for (uint i = 0; i < block.size(); ++i) {
          std::rotate (
            early_mtx[i].begin(), early_mtx[i].begin() + 1, early_mtx[i].end());
        }
        for (uint i = 0; i < block.size(); ++i) {
          early[i][0] += early_mtx[i][0];
          early[i][1] += early_mtx[i][1];
        }
      }
      // late ------------------------------------------------------------------
      for (uint i = 0; i < block.size(); ++i) {
        // the parts between the single sample feedback don't run blockwise
        // unfortunately.
        auto late_mtx_arr = _late_feedback;
        auto late_mtx     = make_crange (late_mtx_arr);

        late_mtx[0] += pre_dif[i][0] * _in_2_late;
        late_mtx[1] += early[i][0] * _er_2_late;
        late_mtx[14] += early[i][1] * _er_2_late;
        late_mtx[15] += pre_dif[i][1] * _in_2_late;

        // diffusion
        auto l = rotation_matrix<8>::tick<float> (
          late_mtx.get_head (8), _late_l_angle);
        crange_copy<float> (late_mtx.get_head (8), l);

        auto r = rotation_matrix<8>::tick<float> (
          late_mtx.advanced (8), _late_r_angle);
        crange_copy<float> (late_mtx.advanced (8), r);

        auto midch = rotation_matrix<8>::tick<float> (
          late_mtx.advanced (4).get_head (8), _late_lr_angle);
        crange_copy<float> (late_mtx.advanced (4).get_head (8), midch);

        late[i][0] = late_mtx[11];
        late[i][1] = late_mtx[4];

        std::rotate (
          late_mtx.begin() + 5, late_mtx.begin() + 6, late_mtx.begin() + 11);

        for (uint i = 0; i < 16; ++i) {
          _late[i].push (make_crange (late_mtx[i]).cast<float_x1>());
        }
        // chorus
        vec<float, 16> n_spls = _late_lfo.tick_filt_sample_and_hold();
        n_spls *= _mod_depth_spls;
        n_spls += _late_n_spls;

        for (uint i = 0; i < 16; ++i) {
          _late_feedback[i]
            = _late[i].get<catmull_rom_interp> (n_spls[i], 0)[0];
        }

        // internal diffusor
        vec<float, 2> mod = _int_dif_lfo.tick_filt_sample_and_hold();
        vec<float, 2> g   = vec_set<2> (_cfg.int_dif.g_base);
        g += mod * _cfg.int_dif.g_mod_depth;
        std::array<float, 2> g_arr {g[0], g[1]};

        std::array<float, 2> diffused;
        diffused[0] = _late_feedback[_cfg.int_dif.channel_l];
        diffused[1] = _late_feedback[_cfg.int_dif.channel_r];

        for (uint st = 0; st < _int_dif.size(); ++st) {
          allpass_stage_tick (
            diffused, _int_dif[st], g_arr, _cfg.int_dif.n_samples[st]);
        }
        _late_feedback[_cfg.int_dif.channel_l] = diffused[0];
        _late_feedback[_cfg.int_dif.channel_r] = diffused[1];

        // Some spice
        auto late_fb_vec
          = vec_tanh_approx_vaneev (vec_from_array (_late_feedback));

        // DC blocker and HP
        late_fb_vec = _filters.tick<dc_idx> (late_fb_vec);

        // Filtering
        auto lp = _filters.tick<lp_idx> (late_fb_vec);
        auto hp = (late_fb_vec - lp) * _filter_hp_att;

        // Final attenuation
        lp *= _rt60_att_l;
        hp *= _rt60_att_h;

        _late_feedback = vec_to_array (lp + hp);
      }
      // output diffusion
      // -----------------------------------------------------------------------
      std::array<float, 2> g_arr {_out_dif_g, _out_dif_g};
      for (uint st = 0; st < _out_dif.size(); ++st) {
        for (uint i = 0; i < block.size(); ++i) {
          allpass_stage_tick (
            late[i], _out_dif[st], g_arr, _cfg.out_dif.n_samples[st]);
        }
      }
      // mixing
      // -----------------------------------------------------------------------
      for (uint i = 0; i < block.size(); ++i) {
        block[i][0] = early[i][0] * _early_gain;
        block[i][1] = early[i][1] * _early_gain;
      }
      for (uint i = 0; i < block.size(); ++i) {
        block[i][0] += late[i][0] * _late_gain;
        block[i][1] += late[i][1] * _late_gain;
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
  void setup_in_diffusor()
  {
    _pre_dif_lfo.reset();
    auto f = vec_set<2> (_cfg.pre_dif.g_mod_freq);
    f[0] *= 1.01; // desync
    _pre_dif_lfo.set_freq (f, _cfg.src.srate);
  }
  //----------------------------------------------------------------------------
  void setup_internal_diffusor()
  {
    _int_dif_lfo.reset();
    auto f = vec_set<2> (_cfg.int_dif.g_mod_freq);
    f[0] *= 1.01; // desync
    _int_dif_lfo.set_freq (f, _cfg.src.srate);
  }
  //----------------------------------------------------------------------------
  void setup_early()
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

    std::reverse (_early_delay_spls[1].begin(), _early_delay_spls[1].end());
    std::reverse (_early_delay_spls[3].begin(), _early_delay_spls[3].end());
  };
  //----------------------------------------------------------------------------
  void setup_late()
  {
    crange_memset<float> (_late_feedback, 0);

    _late_lfo.reset();

    using phase_type = decltype (_late_lfo)::phase_type;
    using value_type = decltype (_late_lfo)::value_type;

    value_type phase;
    float      inc = 1.f / (float) (_late_lfo.n_channels / 2);
    float      ph  = 0.f;
    for (uint i = 0; i < _late_lfo.n_channels / 2; ++i) {
      phase[i]                            = ph;
      phase[i + _late_lfo.n_channels / 2] = ph;
      ph += inc;
    }
    _late_lfo.set_phase (phase_type {phase, phase_type::normalized {}});

    for (uint i = 0; i < _cfg.late.n_samples.size(); ++i) {
      _late_n_spls[i] = (float) (_cfg.late.n_samples[i]);
    }
  }
  //----------------------------------------------------------------------------
  void reset_late_lfo()
  {
    using freq_type             = typename decltype (_late_lfo)::value_type;
    constexpr uint n_side_chnls = vec_traits_t<freq_type>::size / 2;

    auto      freq_fact = _mod_stereo * _cfg.late.max_chorus_width;
    auto      freq_l    = _mod_freq_hz;
    auto      freq_r    = _mod_freq_hz * expf (freq_fact);
    freq_type freq;

    // TODO: "desync_factor": probably a CFG param?
    static constexpr float desync_factor = 0.000001f;
    float                  desync        = 1.f;

    for (uint i = 0; i < n_side_chnls; ++i) {
      freq[i]                = freq_l * desync;
      freq[i + n_side_chnls] = freq_r * desync;
      desync += desync_factor;
    }
    _late_lfo.set_freq (freq, _cfg.src.srate);
  }
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
  template <class T, size_t Sz>
  static uint round_array_pow2_and_accumulate (std::array<T, Sz>& arr, T mod)
  {
    uint mem_total = 0;

    for (uint i = 0; i < arr.size(); ++i) {
      auto spls = arr[i];
      // check constraints for block processing, notice that this is not
      // accounting interpolation
      assert ((spls - mod) > blocksize);
      auto new_spls = pow2_round_ceil (spls + mod);
      assert (new_spls >= spls); // using uint16_t sometimes...
      arr[i] = new_spls;
      mem_total += new_spls;
    }
    return mem_total;
  }
  //----------------------------------------------------------------------------
  void setup_memory()
  {
    // A single contiguous allocation (not likely to matter a lot)
    auto pre_dif_n_spls = _cfg.pre_dif.n_samples;
    auto early_n_spls   = _early_delay_spls;
    auto late_n_spls    = _cfg.late.n_samples;
    auto int_dif_n_spls = _cfg.int_dif.n_samples;
    auto out_dif_n_spls = _cfg.out_dif.n_samples;

    // computing.
    uint mem_total = 0;
    mem_total += round_array_pow2_and_accumulate (pre_dif_n_spls);
    mem_total += round_array_pow2_and_accumulate (early_n_spls);
    mem_total += round_array_pow2_and_accumulate (
      late_n_spls,
      (u16) (_cfg.late.max_chorus_depth_spls + catmull_rom_interp::n_points));
    mem_total += round_array_pow2_and_accumulate (int_dif_n_spls);
    mem_total += round_array_pow2_and_accumulate (out_dif_n_spls);

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
    // late
    for (uint i = 0; i < late_n_spls.size(); ++i) {
      _late[i].reset (mem.cut_head (late_n_spls[i]), 1);
    }
    // internal diffusor
    for (uint i = 0; i < int_dif_n_spls.size(); ++i) {
      for (uint j = 0; j < int_dif_n_spls[0].size(); ++j) {
        _int_dif[i][j].reset (mem.cut_head (int_dif_n_spls[i][j]));
      }
    }
    // output diffusor
    for (uint i = 0; i < out_dif_n_spls.size(); ++i) {
      for (uint j = 0; j < out_dif_n_spls[0].size(); ++j) {
        _out_dif[i][j].reset (mem.cut_head (out_dif_n_spls[i][j]));
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
  static array2d<float, 2, 2> set_rotation_angle (float weight)
  {
    weight *= 0.7f;
    weight += 0.15f; // 0.15 to 0.85

    float w1a = std::cos (weight * 0.5f * M_PI);
    float w1b = std::sqrt (1.f - w1a * w1a);
    float w2a = std::cos (weight * 0.49f * M_PI);
    float w2b = std::sqrt (1.f - w2a * w2a);

    array2d<float, 2, 2> ret;
    ret[0][0] = w1a;
    ret[0][1] = w1b;
    ret[1][0] = w2a;
    ret[1][1] = w2b;
    return ret;
  }
  //----------------------------------------------------------------------------
  void reset_times()
  {
    for (uint i = 0; i < late_cfg::n_channels; ++i) {
      _rt60_att_h[i] = delay_get_feedback_gain_for_time (
        _seconds,
        -60.f,
        _cfg.src.srate,
        vec_set<1> ((float) _cfg.late.n_samples[i]))[0];

      _rt60_att_l[i] = delay_get_feedback_gain_for_time (
        _seconds * _lf_rt60_factor,
        -60.f,
        _cfg.src.srate,
        vec_set<1> ((float) _cfg.late.n_samples[i]))[0];
    }
  }
  //----------------------------------------------------------------------------
  void reset_damping_filters()
  {
    vec<float, 16> freqs = vec_cast<float> (vec_from_array (_cfg.filter.freqs));
    freqs *= (float) std::exp (_damp_freq * _cfg.filter.freq_factor);
    float nyquist = (float) _cfg.src.srate * 0.5f - 2.f;
    freqs         = 1.f - freqs / nyquist; // now a factor ready to be scaled
    // TODO: Why both the damping factor and the damping frequency affect the
    // frequency? this was so on the JSFX prototype but I forgot...
    freqs = nyquist - nyquist * (float) sqrt (_damp_factor) * freqs;
    _filters.reset_coeffs<lp_idx> (freqs, (float) _cfg.src.srate);
    _filter_hp_att = db_to_gain (_cfg.filter.max_att_db * _damp_factor * 0.5f);
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

  std::array<float, late_cfg::n_channels> _late_feedback;
  lfo<late_cfg::n_channels>               _late_lfo;
  array2d<float, 2, 2>                    _late_l_angle;
  array2d<float, 2, 2>                    _late_r_angle;
  array2d<float, 2, 2>                    _late_lr_angle;
  vec<float, late_cfg::n_channels>        _late_n_spls;
  std::array<interpolated_delay_line<float_x1, false>, late_cfg::n_channels>
    _late;

  float _in_2_late      = 1.f;
  float _er_2_late      = 1.f;
  float _mod_freq_hz    = 0.f;
  float _mod_depth_spls = 30.f;
  float _mod_stereo     = 0.f;

  reverb_array<allpass<float_x1>, internal_dif_cfg> _int_dif;
  lfo<2>                                            _int_dif_lfo;

  reverb_array<allpass<float_x1>, internal_dif_cfg> _out_dif;
  float                                             _out_dif_g;

  enum { dc_idx, lp_idx };

  part_classes<
    mp_list<mystran_dc_blocker, onepole_lowpass>,
    vec<float, 16>,
    false>
    _filters;

  float _damp_freq;
  float _damp_factor;
  float _lf_rt60_factor;
  float _filter_hp_att;

  vec<float, late_cfg::n_channels> _rt60_att_h {};
  vec<float, late_cfg::n_channels> _rt60_att_l {};

  float _early_gain = 0.05f;
  float _late_gain  = 0.1f;

  cfg _cfg;

  std::vector<float_x1> _mem;

  float _seconds;
};
//------------------------------------------------------------------------------
} // namespace artv