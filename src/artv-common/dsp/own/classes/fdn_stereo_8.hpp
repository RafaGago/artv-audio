#pragma once

#include "artv-common/misc/compiler.hpp"
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
    float                              g_max;
    array2d<u16, n_channels, n_stages> n_samples;
  };
  //----------------------------------------------------------------------------
  struct early_cfg {
    static constexpr uint n_channels = 4;
    static constexpr uint n_stages   = 4;

    u16   prime_idx;
    float rounding_factor;
    float size_factor;

    struct {
      float meters;
      float span;
      float g;
    } stage[n_stages];
  };
  //----------------------------------------------------------------------------
  struct late_cfg {
    static constexpr uint n_channels = 16;
    static constexpr uint n_stages   = 1;
    float                 max_chorus_width;
    float                 min_chorus_freq;
    float                 max_chorus_freq;
    // max mod depth point (and below), linearly decreasing after that
    float max_chorus_depth_freq;
    float size_factor; // spls = [-1,1] * exp (size_factor)
    u16   max_chorus_depth_spls;
    u16   prime_idx;
    std::array<u16, n_channels> n_samples;
  };
  //----------------------------------------------------------------------------
  struct internal_dif_cfg {
    static constexpr uint              n_channels = 2;
    static constexpr uint              n_stages   = 4;
    float                              g_mod_depth;
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
  struct stereo_cfg {
    static constexpr uint n_channels = 2;
    float                 freq_center;
    float                 freq_factor;
    float                 g_base;
    float                 max_samples;
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
    stereo_cfg       stereo;
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
    r.pre_dif.g_max       = 0.39f;

    r.pre_dif.n_samples = make_array (
      array_cast<u16> (make_array (43, 43)),
      array_cast<u16> (make_array (113, 113)),
      array_cast<u16> (make_array (317, 317)),
      array_cast<u16> (make_array (907, 906)));

    r.early.stage[0].meters = 4.129f;
    r.early.stage[0].span   = 3.209f;
    r.early.stage[0].g      = 0.5f;
    r.early.stage[1].meters = 5.101f;
    r.early.stage[1].span   = 5.003f;
    r.early.stage[1].g      = 0.3f;
    r.early.stage[2].meters = 3.163f;
    r.early.stage[2].span   = 3.121f;
    r.early.stage[2].g      = 0.5f;
    r.early.stage[3].meters = 7.333f;
    r.early.stage[3].span   = 3.187f;
    r.early.stage[3].g      = 0.65f;
    r.early.prime_idx       = 1;
    r.early.rounding_factor = 1000;
    r.early.size_factor     = 1.f; // log(e)

    r.late.prime_idx   = 15;
    r.late.size_factor = 2.2f;

    r.late.max_chorus_freq       = 4.5f;
    r.late.min_chorus_freq       = 0.15f;
    r.late.max_chorus_depth_spls = 150; // bipolar, 2x the samples here
    r.late.max_chorus_depth_freq = 0.15f;
    r.late.max_chorus_width      = 0.15f;

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

    r.filter.max_att_db  = -13.f;
    r.filter.freq_factor = std::log (4.5f);

    r.filter.freqs = array_cast<u16> (make_array (
      3000,
      2800,
      2000,
      2200,
      1000,
      1200,
      1100,
      1000,
      900,
      940,
      800,
      860,
      600,
      640,
      330,
      300));

    from_ascending_pairs_to_internal_chnl_order (r.filter.freqs);

    r.stereo.freq_factor = 0.317f;
    r.stereo.g_base      = 0.4f;
    r.stereo.max_samples = 20.f;

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
  void set_early_gain (float g)
  {
    static constexpr float zero = constexpr_db_to_gain (-60.f);
    _early_gain                 = g > zero ? g : 0.f;
  }
  //----------------------------------------------------------------------------
  void set_late_gain (float g)
  {
    static constexpr float zero = constexpr_db_to_gain (-60.f);
    _late_gain                  = g > zero ? g : 0.f;
  }
  //----------------------------------------------------------------------------
  void set_mod_freq (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    auto diff    = _cfg.late.max_chorus_freq - _cfg.late.min_chorus_freq;
    _mod_freq_hz = _cfg.late.min_chorus_freq + (diff * factor * factor);
    reset_mod_freq();
    reset_mod_depth();
  }
  //----------------------------------------------------------------------------
  void set_mod_depth (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _mod_depth_factor = factor * factor;
    reset_mod_depth();
  }
  //----------------------------------------------------------------------------
  void set_mod_stereo (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    _mod_stereo = factor;
    reset_mod_freq();
    reset_mod_stereo();
  }
  //----------------------------------------------------------------------------
  void set_mod_wave (uint wv)
  {
    assert (wv < modwv_count);
    _late_wave = wv;
    reset_mod_freq();
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
    assert (factor >= -1.f && factor <= 1.f);
    _late_l_angle = set_rotation_angle (factor);
  }
  //----------------------------------------------------------------------------
  void set_r_matrix_angle (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    _late_r_angle = set_rotation_angle (factor);
  }
  //----------------------------------------------------------------------------
  void set_lr_matrix_angle (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    _late_lr_angle = set_rotation_angle (factor);
  }
  //----------------------------------------------------------------------------
  void set_time_msec (float msec)
  {
    _seconds = msec * 0.001f;
    reset_times();
  }
  //----------------------------------------------------------------------------
  void set_size (float size)
  {
    assert (size >= -1.f && size <= 1.f);
    _size = size;
    reset_times();
    reset_mod_depth();
  }
  //----------------------------------------------------------------------------
  void set_damp_freq (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    vec<float, 16> freqs = vec_cast<float> (vec_from_array (_cfg.filter.freqs));
    freqs *= (float) std::exp (factor * _cfg.filter.freq_factor);
    _filters.reset_coeffs<lp_idx> (freqs, (float) _cfg.src.srate);
  }
  //----------------------------------------------------------------------------
  void set_damp_factor (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    factor *= factor;
    _filter_hp_att = db_to_gain (_cfg.filter.max_att_db * factor * 0.5f);
  }
  //----------------------------------------------------------------------------
  void set_hp_freq (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    float hp_freq = exp (factor * 4.1f); // TODO: make variables
    _filters.reset_coeffs<dc_idx> (
      vec_set<16> (hp_freq), (float) _cfg.src.srate);
  }
  //----------------------------------------------------------------------------
  void set_lf_time_factor (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    _lf_rt60_factor = exp (factor); // TODO: make variables
    reset_times();
  }
  //----------------------------------------------------------------------------
  void set_predelay (float sixteenths)
  {
    assert (sixteenths >= 0.f && sixteenths <= 16.f);
    float predelay_spls = sixteenths * _beat_16th_spls;
    predelay_spls       = predelay_spls - (float) _pre_delay_lat_spls;
    _pre_delay_spls     = predelay_spls <= 0.f ? 0 : (uint) predelay_spls;
  }
  //----------------------------------------------------------------------------
  void set_gap (float sixteenths)
  {
    assert (sixteenths >= 0.f && sixteenths <= 16.f);
    _gap_spls = (uint) sixteenths * _beat_16th_spls;
  }
  //----------------------------------------------------------------------------
  void set_er_size (float factor)
  {
    assert (factor >= -1.f && factor <= 1.f);
    reset_early (factor);
  }
  //----------------------------------------------------------------------------
  void set_stereo (float factor)
  {
    assert (factor >= 0.f && factor <= 1.f);
    _stereo = sqrt (factor);
  }
  //----------------------------------------------------------------------------
  void set_test_param (uint v) { _test = v; }
  //----------------------------------------------------------------------------
  void reset (uint samplerate, float bpm, cfg const& cfg)
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

    _pre_delay_lat_spls = get_in_diff_delay();
    float sixtenths_sec = bpm * (1.f / 60.f) * (1.f / 60.f);
    _beat_16th_spls     = sixtenths_sec * (float) samplerate;

    setup_late();
    setup_memory();
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
  uint get_in_diff_delay() const
  {
    uint ret = 0;
    for (auto v : _cfg.pre_dif.n_samples) {
      ret += std::max (v[0], v[0]);
    }
    return ret;
  }
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

    static_assert (early_cfg::n_channels == 4); // hardcoding for readability
    array2d<float, 4, blocksize> early_mtx;

    static_assert (late_cfg::n_channels == 16); // hardcoding for readability
    array2d<float, 2, blocksize>                           early_then_late;
    alignas (vec<float, 16>) array2d<float, 16, blocksize> late_mtx;

    while (io.size()) {
      auto block = io.cut_head (std::min<uint> (io.size(), blocksize));

      // pre diffusor ----------------------------------------------------------

      // AP gain LFO
      auto mod_g     = make_crange (tmp);
      uint late_wave = _late_wave; // telling the optimizer to ignore changes
      for (uint i = 0; i < block.size(); ++i) {
        vec<float, 2> mod;
        if (late_wave == modwv_sh) {
          mod = _int_dif_lfo.tick_filt_sample_and_hold();
        }
        else {
          mod = _int_dif_lfo.tick_sine();
        }
        vec<float, 2> g = vec_set<2> (_pre_dif_g);
        g += mod * _cfg.pre_dif.g_mod_depth;
        mod_g[i] = vec_to_array (g);
      }
      // hardcoding 2's to avoid unreadability, so "static_asserting"
      static_assert (n_channels == 2);
      if (_pre_delay_spls > _pre_delay_lat_spls) {
        uint pre_delay_spls = _pre_delay_spls - _pre_delay_lat_spls;
        for (uint i = 0; i < block.size(); ++i) {
          _pre_delay.push (block[i]);
          pre_dif[i][0] = _pre_delay.get (pre_delay_spls, 0);
          pre_dif[i][1] = _pre_delay.get (pre_delay_spls, 1);
        }
      }
      else {
        crange_copy<std::array<float, 2>> (pre_dif, block);
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

      // building the 4-wide matrices
      for (uint i = 0; i < block.size(); ++i) {
        early_mtx[i][0] = pre_dif[i][0];
        early_mtx[i][1] = pre_dif[i][1];
        early_mtx[i][2] = (early_mtx[i][0] + early_mtx[i][1]) * 0.5f; // mid
        early_mtx[i][3] = early_mtx[i][2];

        block[i][0] = 0.f;
        block[i][1] = 0.f;
      }

      if (_early_gain != 0.f || _er_2_late != 0.f) {
        float stage_gain = 0.7f;
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
              early_mtx[i] = hadamard_matrix<4>::tick<float> (early_mtx[i]);
            }
          }
          else {
            for (uint i = 0; i < block.size(); ++i) {
              early_mtx[i] = householder_matrix<4>::tick<float> (early_mtx[i]);
            }
          }
          // rotation + storage
          for (uint i = 0; i < block.size(); ++i) {
            std::rotate (
              early_mtx[i].begin(),
              early_mtx[i].begin() + 1,
              early_mtx[i].end());
            block[i][(stage & 1) == 0] += early_mtx[i][0] * stage_gain;
            block[i][(stage & 1) == 1] += early_mtx[i][1] * stage_gain;
          }
          stage_gain *= stage_gain;
        }
      }
      else {
        crange_memset (block, 0);
      }
      // gap -------------------------------------------------------------------
      if (_gap_spls > _gap_lat_spls) {
        uint gap_spls = _gap_spls - _gap_lat_spls;

        for (uint i = 0; i < block.size(); ++i) {
          // reminder block[i][0/1] contains the unscaled early reflections
          std::array<float, 4> gap_spl {
            pre_dif[i][0], pre_dif[i][1], block[i][0], block[i][1]};

          _gap.push (gap_spl);

          pre_dif[i][0]         = _gap.get (gap_spls, 0);
          pre_dif[i][1]         = _gap.get (gap_spls, 1);
          early_then_late[i][0] = _gap.get (gap_spls, 2);
          early_then_late[i][1] = _gap.get (gap_spls, 3);
        }
      }
      else {
        crange_copy<std::array<float, 2>> (early_then_late, block);
      }
      // late
      // -----------------------------------------------------------------------
      if (_late_gain != 0.f) {

        // feedback chorus lfo
        uint late_wave = _late_wave; // telling the optimizer to ignore changes
        for (uint i = 0; i < block.size(); ++i) {
          vec<float, 16> mod;
          switch (late_wave) {
          case modwv_sh:
            mod = _late_lfo.tick_filt_sample_and_hold();
            break;
          case modwv_sin:
            mod = _late_lfo.tick_sine();
            break;
          case modwv_tri:
            mod = _late_lfo.tick_triangle();
            break;
          case modwv_tra:
            mod = _late_lfo.tick_trapezoid (vec_set<16> (0.75f));
            break;
          default:
            assert (false);
            break;
          };
          auto n_spls_mod = mod * _mod_depth_spls;
          auto n_spls     = n_spls_mod + _late_n_spls;
          for (uint j = 0; j < 16; ++j) {
            assert (n_spls[j] >= 0.f);
          }
          late_mtx[i] = _late.get (vec_to_array (n_spls - (float) i));
        }
        // internal diffusor lfo.
        mod_g = make_crange (tmp);
        for (uint i = 0; i < block.size(); ++i) {
          vec<float, 2> mod;
          if (late_wave == modwv_sh) {
            mod = _int_dif_lfo.tick_filt_sample_and_hold();
          }
          else {
            mod = _int_dif_lfo.tick_sine();
          }
          vec<float, 2> g = vec_set<2> (_cfg.int_dif.g_base);
          g += mod * _cfg.int_dif.g_mod_depth;
          mod_g[i] = vec_to_array (g);
        }

        // internal diffusor
        for (uint i = 0; i < block.size(); ++i) {
          std::array<float, 2> diffused;
          diffused[0] = late_mtx[i][_cfg.int_dif.channel_l];
          diffused[1] = late_mtx[i][_cfg.int_dif.channel_r];

          for (uint st = 0; st < _int_dif.size(); ++st) {
            allpass_stage_tick (
              diffused, _int_dif[st], mod_g[i], _cfg.int_dif.n_samples[st]);
          }
          late_mtx[i][_cfg.int_dif.channel_l] = diffused[0];
          late_mtx[i][_cfg.int_dif.channel_r] = diffused[1];
        }

        // sigmoid + dc + filtering in one stage to avoid multiple vector to
        // SIMD conversions
        for (uint i = 0; i < block.size(); ++i) {
          auto late_fb_vec = vec_from_array (late_mtx[i]);
          late_fb_vec
            = late_fb_vec / vec_sqrt (late_fb_vec * late_fb_vec + 1.f);
          late_fb_vec = _filters.tick<dc_idx> (late_fb_vec);
          auto lp     = _filters.tick<lp_idx> (late_fb_vec);
          auto hp     = (late_fb_vec - lp) * _filter_hp_att;
          // Final attenuation
          lp *= _rt60_att_l;
          hp *= _rt60_att_h;
          late_mtx[i] = vec_to_array (lp + hp);
        }

        // end of feedback path. Starting feedforward
        for (uint i = 0; i < block.size(); ++i) {
          // reminder "early_then_late" contains the unscaled early reflections
          late_mtx[i][0] += pre_dif[i][0] * _in_2_late;
          late_mtx[i][1] += early_then_late[i][1] * _er_2_late;
          late_mtx[i][14] += early_then_late[i][0] * _er_2_late;
          late_mtx[i][15] += pre_dif[i][1] * _in_2_late;
        }

        for (uint i = 0; i < block.size(); ++i) {
          auto lm = make_crange (late_mtx[i]);
          // diffusion
          auto l
            = rotation_matrix<8>::tick<float> (lm.get_head (8), _late_l_angle);
          crange_copy<float> (lm.get_head (8), l);

          auto r
            = rotation_matrix<8>::tick<float> (lm.advanced (8), _late_r_angle);
          crange_copy<float> (lm.advanced (8), r);

          auto midch = rotation_matrix<8>::tick<float> (
            lm.advanced (4).get_head (8), _late_lr_angle);
          crange_copy<float> (lm.advanced (4).get_head (8), midch);
        }

        for (uint i = 0; i < block.size(); ++i) {
          early_then_late[i][0] = late_mtx[i][11];
          early_then_late[i][1] = late_mtx[i][4];

          std::rotate (
            late_mtx[i].begin() + 5,
            late_mtx[i].begin() + 6,
            late_mtx[i].begin() + 11);

          _late.push (late_mtx[i]);
        }
        // stereo comb
        // ---------------------------------------------------------------------
        uint st_delay
          = _cfg.stereo.max_samples * _mod_depth_factor * _mod_depth_factor;
        st_delay += blocksize;
        for (uint i = 0; i < block.size(); ++i) {
          vec<float, 1> stmod;
          if (late_wave == modwv_sh) {
            stmod = _stereo_lfo.tick_filt_sample_and_hold();
          }
          else {
            stmod = _stereo_lfo.tick_sine();
          }
          stmod *= _cfg.stereo.g_base * _mod_stereo;

          early_then_late[i][0] = _stereo_allpass[0].tick (
            make_vec (early_then_late[i][0]),
            st_delay,
            vec_set<1> (0.f),
            stmod)[0];
          early_then_late[i][1] = _stereo_allpass[1].tick (
            make_vec (early_then_late[i][1]),
            st_delay,
            vec_set<1> (0.f),
            -stmod)[0];
        }
        // output diffusion
        // ---------------------------------------------------------------------
        std::array<float, 2> g_arr {_out_dif_g, _out_dif_g};
        for (uint st = 0; st < _out_dif.size(); ++st) {
          for (uint i = 0; i < block.size(); ++i) {
            allpass_stage_tick (
              early_then_late[i],
              _out_dif[st],
              g_arr,
              _cfg.out_dif.n_samples[st]);
          }
        }
        // mixing late + early + stereo
        for (uint i = 0; i < block.size(); ++i) {
          // previously early had no gain scaling (optimization to remove the
          // need for an intermediate buffer)
          block[i][0] *= _early_gain;
          block[i][1] *= _early_gain;

          block[i][0] += early_then_late[i][0] * _late_gain;
          block[i][1] += early_then_late[i][1] * _late_gain;

          block[i][1] = block[i][1] * _stereo + block[i][0] * (1.f - _stereo);
        }
      }
      else {
        // scale the early reflections + stereo
        for (uint i = 0; i < block.size(); ++i) {
          block[i][0] *= _early_gain;
          block[i][1] *= _early_gain;

          block[i][1] = block[i][1] * _stereo + block[i][0] * (1.f - _stereo);
        }
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
  void reset_early (float size_factor)
  {
    float f             = exp (size_factor * _cfg.early.size_factor);
    _pre_delay_lat_spls = (uint) -1;
    // compute delay lengths
    for (uint i = 0; i < _cfg.early.n_stages; ++i) {
      uint spls_min = delay_length::meters_to_samples (
        _cfg.early.stage[i].meters * f, _cfg.src.srate, false);
      uint spls_max = delay_length::meters_to_samples (
        _cfg.early.stage[i].meters * f * _cfg.early.stage[i].span,
        _cfg.src.srate,
        true);
      get_delay_length<u16> (
        _early_delay_spls[i],
        spls_min,
        spls_max,
        _cfg.early.prime_idx,
        _cfg.early.rounding_factor);
      _pre_delay_lat_spls = std::min<uint> (_pre_delay_lat_spls, spls_min);
    }

    std::reverse (_early_delay_spls[1].begin(), _early_delay_spls[1].end());
    std::reverse (_early_delay_spls[3].begin(), _early_delay_spls[3].end());
  };
  //----------------------------------------------------------------------------
  void setup_late()
  {
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
      _late_n_spls_master[i] = (float) (_cfg.late.n_samples[i]);
    }
#if 0
    for (auto& dl : _late) {
      dl.set_resync_delta (20.f);
    }
#endif
  }
  //----------------------------------------------------------------------------
  void reset_mod_freq()
  {
    using vec_type              = typename decltype (_late_lfo)::value_type;
    constexpr uint n_side_chnls = vec_traits_t<vec_type>::size / 2;
    constexpr uint n_chnls      = vec_traits_t<vec_type>::size;

    vec_type vfreq {};

    auto freq_fact = _mod_stereo * _cfg.late.max_chorus_width;
    auto freq_l    = _mod_freq_hz;
    auto freq_r    = _mod_freq_hz * expf (freq_fact);

    if (_late_wave == modwv_sh) {
      float desync_factor = 0.08f * _mod_stereo;
      float desync        = 1.f;

      for (uint i = 0; i < n_side_chnls; ++i) {
        vfreq[i] = freq_l * desync;
        desync += desync_factor;
      }
      desync = 1.f;
      for (uint i = 0; i < n_side_chnls; ++i) {
        vfreq[i + n_side_chnls] = freq_r * desync;
        desync += desync_factor;
      }
    }
    else {
      for (uint i = 0; i < n_side_chnls; ++i) {
        vfreq[i] = freq_l;
      }
      for (uint i = 0; i < n_side_chnls; ++i) {
        vfreq[i + n_side_chnls] = freq_r;
      }
    }

    _late_lfo.set_freq (vfreq, _cfg.src.srate);

    vec<float, 2> dif_freq {freq_r, freq_l};
    _int_dif_lfo.set_freq (dif_freq, _cfg.src.srate);

    _pre_dif_lfo.set_freq (dif_freq, _cfg.src.srate);

    _stereo_lfo.set_freq (
      make_vec (freq_r + freq_l * _cfg.stereo.freq_factor), _cfg.src.srate);
  }
  //----------------------------------------------------------------------------
  void reset_mod_stereo()
  {
    using vec_type              = typename decltype (_late_lfo)::value_type;
    constexpr uint n_side_chnls = vec_traits_t<vec_type>::size / 2;
    constexpr uint n_chnls      = vec_traits_t<vec_type>::size;

    vec_type vphase {};

    if (_late_wave == modwv_sh) {
      float inc    = 1.f / ((float) (n_chnls));
      float cphase = 0.f;
      inc *= _mod_stereo;
      for (uint i = 0; i < n_chnls; ++i) {
        vphase[i] = cphase;
        cphase += inc;
      }
    }
    else {
      float ph = 0.5 * _mod_stereo;

      for (uint i = 0; i < n_side_chnls; ++i) {
        vphase[i] = (i & 2) ? ph : -ph;
      }
      for (uint i = 0; i < n_side_chnls; ++i) {
        vphase[i] = (i & 2) ? -ph : ph;
      }
    }
    _late_lfo.set_phase (phase<16> {vphase, phase<16>::normalized()});
  }
  //----------------------------------------------------------------------------
  template <class T, size_t A, size_t B>
  static uint convert_to_max_sizes (array2d<T, A, B>& arr)
  {
    uint mem_total = 0;
    for (uint i = 0; i < B; ++i) {
      for (uint j = 0; j < A; ++j) {
        auto spls = arr[i][j];
        // check constraints for block processing, notice that this is not
        // accounting interpolation
        assert (spls > blocksize);
        auto new_spls = pow2_round_ceil (spls + 1);
        assert (new_spls >= spls); // using uint16_t sometimes...
        arr[i][j] = new_spls;
        mem_total += new_spls;
      }
    }
    return mem_total;
  }
  //----------------------------------------------------------------------------
  template <class T, size_t Sz>
  static uint convert_to_max_sizes (std::array<T, Sz>& arr, T mod_add)
  {
    uint mem_total = 0;

    for (uint i = 0; i < arr.size(); ++i) {
      auto spls = arr[i];
      // check constraints for block processing, notice that this is not
      // accounting interpolation
      assert ((spls - mod_add) > blocksize);
      auto new_spls = pow2_round_ceil (spls + mod_add + 1);
      assert (new_spls >= spls); // using uint16_t sometimes...
      arr[i] = new_spls;
      mem_total += new_spls;
    }
    return mem_total;
  }
  //----------------------------------------------------------------------------
  void setup_memory()
  {
    // set "_early_delay_spls" to its maximum possible size just for this
    // maximum memory computation.
    reset_early (1.f);

    // A single contiguous allocation (not likely to matter a lot)
    auto pre_dif_sizes = _cfg.pre_dif.n_samples;
    auto early_sizes   = _early_delay_spls;
    auto int_dif_sizes = _cfg.int_dif.n_samples;
    auto out_dif_sizes = _cfg.out_dif.n_samples;

    auto late_sizes_tmp = _late_n_spls_master;
    for (auto& spls : late_sizes_tmp) {
      spls *= exp (_cfg.late.size_factor);
      spls += 1.f; // before casting
    }
    auto late_sizes = array_cast<u32> (late_sizes_tmp);

    // computing.
    uint mem_total = 0;
    mem_total += convert_to_max_sizes (pre_dif_sizes);
    mem_total += convert_to_max_sizes (early_sizes);
    mem_total += convert_to_max_sizes (
      late_sizes,
      (u32) (_cfg.late.max_chorus_depth_spls + delay_line_type::interp::n_points));
#if 0
    for (uint i = 0; i < late_sizes.size(); ++i) {
      late_sizes[i] += _late[i].interp_overhead_elems (1);
    }
    mem_total += _late[0].interp_overhead_elems (1) * late_sizes.size();
#endif

    mem_total += convert_to_max_sizes (int_dif_sizes);
    mem_total += convert_to_max_sizes (out_dif_sizes);

    // pre_delay
    uint pre_delay_spls_max = (uint) (_beat_16th_spls * 16.f + 1.f);
    pre_delay_spls_max -= _pre_delay_lat_spls;
    pre_delay_spls_max *= 2;
    pre_delay_spls_max = pow2_round_ceil (pre_delay_spls_max);
    mem_total += pre_delay_spls_max;

    // gap
    uint gap_spls_max = (uint) (_beat_16th_spls * 16.f + 1.f);
    gap_spls_max *= 4;
    gap_spls_max = pow2_round_ceil (gap_spls_max);
    mem_total += gap_spls_max;

    // stereo
    uint stereo_allpass_size
      = pow2_round_ceil ((uint) _cfg.stereo.max_samples + blocksize);
    mem_total += stereo_allpass_size * 2;

    // allocating_mem
    _mem.clear();
    _mem.resize (mem_total);

    // assigning memory
    auto mem = make_crange (_mem);
    // pre delay
    _pre_delay.reset (mem.cut_head (pre_delay_spls_max).cast (float {}), 2);
    // gap
    _gap.reset (mem.cut_head (gap_spls_max).cast (float {}), 4);
    // pre diffusor
    for (uint i = 0; i < pre_dif_sizes.size(); ++i) {
      for (uint j = 0; j < pre_dif_sizes[0].size(); ++j) {
        _pre_dif[i][j].reset (mem.cut_head (pre_dif_sizes[i][j]));
      }
    }
    // early
    for (uint i = 0; i < early_sizes.size(); ++i) {
      for (uint j = 0; j < early_sizes[0].size(); ++j) {
        _early[i][j].reset (mem.cut_head (early_sizes[i][j]));
      }
    }
    // late
    mem = _late.reset<u32> (mem.cast<float>(), late_sizes).cast<float_x1>();
    // internal diffusor
    for (uint i = 0; i < int_dif_sizes.size(); ++i) {
      for (uint j = 0; j < int_dif_sizes[0].size(); ++j) {
        _int_dif[i][j].reset (mem.cut_head (int_dif_sizes[i][j]));
      }
    }
    // output diffusor
    for (uint i = 0; i < out_dif_sizes.size(); ++i) {
      for (uint j = 0; j < out_dif_sizes[0].size(); ++j) {
        _out_dif[i][j].reset (mem.cut_head (out_dif_sizes[i][j]));
      }
    }
    // stereo
    for (uint i = 0; i < _stereo_allpass.size(); ++i) {
      _stereo_allpass[i].reset (mem.cut_head (stereo_allpass_size));
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
    // from [-1,1] to [0, 1]
    weight += 1.f;
    weight *= 0.5f;
    // from [0,1] to [0.1, 0.9]
    weight *= 0.8f;
    weight += 0.1f;

    float w1a = cos (weight * 0.5f * M_PI);
    float w1b = sqrt (1.f - w1a * w1a);
    float w2a = cos (weight * 0.49f * M_PI);
    float w2b = sqrt (1.f - w2a * w2a);

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
    _late_n_spls = vec_from_array (_late_n_spls_master) * (float) exp (_size);

    _gap_lat_spls = (uint) -1;

    for (uint i = 0; i < late_cfg::n_channels; ++i) {
      _rt60_att_h[i] = delay_get_feedback_gain_for_time (
        _seconds, -60.f, _cfg.src.srate, vec_set<1> (_late_n_spls[i]))[0];

      _rt60_att_l[i] = delay_get_feedback_gain_for_time (
        _seconds * _lf_rt60_factor,
        -60.f,
        _cfg.src.srate,
        vec_set<1> (_late_n_spls[i]))[0];

      _gap_lat_spls = std::min<uint> (_late_n_spls[i], _gap_lat_spls);
    }
  }
  //----------------------------------------------------------------------------
  void reset_mod_depth()
  {
    // kind-of constant module
    float spls = (float) _cfg.late.max_chorus_depth_spls;
    spls *= _mod_depth_factor * _cfg.late.max_chorus_depth_freq;
    spls /= _mod_freq_hz;

    // limit the excursion on both sides
    auto max_spls = (uint) _late_n_spls[0]
      - (blocksize + delay_line_type::interp::n_points);

    _mod_depth_spls = std::min<float> (spls, max_spls);
    _mod_depth_spls
      = std::min<float> (_mod_depth_spls, _cfg.late.max_chorus_depth_spls);
  }
  //----------------------------------------------------------------------------
  template <class T, class Cfg>
  using reverb_array = array2d<T, Cfg::n_channels, Cfg::n_stages>;
  //----------------------------------------------------------------------------
  enum mod_wave {
    modwv_sh,
    modwv_sin,
    modwv_tri,
    modwv_tra,
    modwv_count,
  };

  block_resampler<float, 2> _resampler;

  static_delay_line<float, true, true> _pre_delay;
  uint                                 _pre_delay_lat_spls;
  uint                                 _pre_delay_spls;

  reverb_array<allpass<float_x1>, pre_dif_cfg> _pre_dif;
  lfo<2>                                       _pre_dif_lfo;
  float                                        _pre_dif_g;

  reverb_array<u16, early_cfg>               _early_delay_spls;
  reverb_array<allpass<float_x1>, early_cfg> _early;

  uint                                 _gap_lat_spls;
  uint                                 _gap_spls;
  static_delay_line<float, true, true> _gap;

  lfo<late_cfg::n_channels>               _late_lfo;
  array2d<float, 2, 2>                    _late_l_angle;
  array2d<float, 2, 2>                    _late_r_angle;
  array2d<float, 2, 2>                    _late_lr_angle;
  std::array<float, late_cfg::n_channels> _late_n_spls_master;
  vec<float, late_cfg::n_channels>        _late_n_spls;
  uint                                    _late_wave = 0;

  using delay_line_type = interpolated_multisized_delay_line<
    float,
    late_cfg::n_channels,
    catmull_rom_interp,
    true>;
  delay_line_type _late;

  float _size             = 1.f;
  float _in_2_late        = 1.f;
  float _er_2_late        = 1.f;
  float _mod_freq_hz      = 0.f;
  float _mod_depth_spls   = 30.f;
  float _mod_depth_factor = 0.f;
  float _mod_stereo       = 0.f;

  reverb_array<allpass<float_x1>, internal_dif_cfg> _int_dif;
  lfo<2>                                            _int_dif_lfo;

  reverb_array<allpass<float_x1>, internal_dif_cfg> _out_dif;
  float                                             _out_dif_g;

  std::array<allpass<float_x1>, n_channels> _stereo_allpass;
  lfo<1>                                    _stereo_lfo;

  enum { dc_idx, lp_idx };

  part_classes<
    mp_list<mystran_dc_blocker, onepole_lowpass>,
    vec<float, 16>,
    false>
    _filters;

  float _lf_rt60_factor;
  float _filter_hp_att;

  vec<float, late_cfg::n_channels> _rt60_att_h {};
  vec<float, late_cfg::n_channels> _rt60_att_l {};

  float _early_gain = 0.05f;
  float _late_gain  = 0.1f;

  cfg _cfg;

  std::vector<float_x1> _mem;

  float _early_size;
  float _seconds;
  float _stereo;
  float _beat_16th_spls;

  uint _test;
};
//------------------------------------------------------------------------------
} // namespace artv
