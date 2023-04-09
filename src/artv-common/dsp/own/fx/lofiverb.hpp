#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/ducker.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/fx/lofiverb-engine.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/own/parts/waveshapers/sigmoid.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

// #define LOFIVERB_ADD_DEBUG_ALGO 1

namespace artv {
namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
static constexpr uint max_block_size = 32;

#ifdef LOFIVERB_ADD_DEBUG_ALGO
//------------------------------------------------------------------------------
static constexpr auto get_debug_algo_spec()
{
  return make_array<stage_data> (
    make_ap (147, 0.5), // 0
    make_ap (183, 0.4), // 1
    make_ap (389, 0.3), // 2
    make_block_delay (max_block_size) // 3
  );
}

struct debug_algo_spec {
  static constexpr auto values {get_debug_algo_spec()};
};
#endif
//------------------------------------------------------------------------------
static constexpr auto get_abyss_spec()
{
  return make_array<stage_data> (
    // diffusors
    make_ap (147, -0.707), // 0
    make_ap (183, 0.707), // 1
    make_ap (389, -0.6), // 2
    make_ap (401, 0.6), // 3
    // er (not ER at the end, as this is now another reverb...)
    make_ap (1367, 0.35, 71 + 70), // 4
    make_lp(), // 5
    make_lp(), // 6
    make_quantizer(), // 7
    make_ap (1787, 0.5, 261), // 8
    make_block_delay (max_block_size), // 9 to allow block processing
    // loop1
    make_ap (977, 0.5 /*overridden*/, 51), // 10
    make_delay (2819), // 11
    make_lp(), // 12
    make_lp(), // 13
    make_quantizer(), // 14
    make_ap (863, -0.5 /*overridden*/), // 15
    make_delay (1021), // 16
    make_quantizer(), // 17
    make_ap (1453, 0.618), // 18
    make_block_delay (
      787), // 19 delay (allows block processing) (> blocksz + 1)
    // loop2
    make_ap (947, 0.5 /*overridden*/, 67), // 20
    make_delay (3191), // 21
    make_lp(), // 22
    make_lp(), // 23
    make_quantizer(), // 24
    make_ap (887, -0.5 /*overridden*/), // 25
    make_delay (1049), // 26
    make_quantizer(), // 27
    make_ap (1367, 0.618), // 28
    make_hp (0.98), // 29
    make_block_delay (
      647)); // 30 delay (allows block processing) (> blocksz + 1)
}

struct abyss_spec {
  static constexpr auto values {get_abyss_spec()};
};
//------------------------------------------------------------------------------
static constexpr auto get_small_space_spec()
{
  return make_array<stage_data> (
    // left channel output conditioning
    make_ap (71), // 0
    make_ap (116), // 1
    make_ap (172), // 2
    make_ap (277), // 3

    // right channel output conditioning
    make_ap (72), // 4
    make_ap (117), // 5
    make_ap (175), // 6
    make_ap (274), // 7

    // m diffusors
    make_ap (23, 0.8), // 8
    make_ap (130, 0.8), // 9
    make_ap (217, 0.8), // 10
    // s diffusors
    make_ap (22, 0.8), // 11
    make_ap (131, 0.8), // 12
    make_ap (219, 0.8), // 13

    make_quantizer(), // 14
    make_quantizer(), // 15
    make_quantizer(), // 16
    make_quantizer(), // 17
    // block a iteration 1
    make_ap (153), // nested 3x (mod g) // 18
    make_ap (89, -0.04), // 19
    make_ap (60, 0.04), // 20
    make_delay (201), // 21

    // block b iteration 1
    make_ap (139), // nested 3x (mod g) // 22
    make_ap (79, -0.04), // 23
    make_ap (53, 0.04), // 24
    make_delay (185), // mod delay // 25

    // block c iteration 1
    make_ap (149), // nested 3x (mod g) // 26
    make_ap (83, 0.04), // 27
    make_ap (59, -0.04), // 28
    make_delay (193), // 29

    // block d iteration 1
    make_ap (167), // nested 3x (mod g) // 30
    make_ap (97, -0.04), // 31
    make_ap (67, 0.04), // 32
    make_delay (221), // mod delay // 33

    make_quantizer(), // 34
    make_quantizer(), // 35
    make_quantizer(), // 36
    make_quantizer(), // 37

    // block a iteration 2
    make_lp(), // 38
    make_lp(), // 39
    make_ap (113, 0.1), // nested 3x // 40
    make_ap (67, 0.04), // 41
    make_ap (47, -0.04), // nested end // 42
    make_delay (122), // 43
    make_ap (119, 0.6), // 44
    make_ap (67, 0.6), // 45
    make_ap (47, -0.6), // 46
    make_block_delay (32), // feedback point // 47

    // block b iteration 2
    make_lp(), // 48
    make_lp(), // 49
    make_ap (114, -0.1), // nested 3x // 50
    make_ap (66, -0.04), // 51
    make_ap (47, -0.04), // 52
    make_ap (9, 0, 2), // single. modulated both gain and time // 53
    make_block_delay (149), // feedback point // 54

    // block c iteration 2
    make_lp(), // 55
    make_lp(), // 56
    make_ap (116, -0.1), // nested 3x // 57
    make_ap (65, -0.04), // 58
    make_ap (46, -0.04), // 59
    make_ap (9, 0, 3), // single. modulated both gain and time // 60
    make_block_delay (151), // feedback point // 61

    // block d iteration 2
    make_lp(), // 62
    make_lp(), // 63
    make_ap (121, -0.1), // nested 3x // 64
    make_ap (69, 0.04), // 65
    make_ap (47, 0.04), // 66
    make_block_delay (157) // feedback point // 67
  );
}

struct small_space_spec {
  static constexpr auto values {get_small_space_spec()};
};
//------------------------------------------------------------------------------
static constexpr auto get_midifex49_spec()
{
  return make_array<stage_data> (
    // diffusors
    make_ap (321, 0.5), // 0 PreAP
    make_ap (431, 0.5), // 1 PreAP
    make_ap (968, 0.5), // 2 PreAP
    make_ap (1620, 0.5), // 3 PreAP

    make_delay (21), // 4 L1
    make_block_delay (1010, 0), // 5 R1

    make_delay (1624), // 6 Loop
    make_ap (1992, 0.5, 17), // 7 Loop

    make_block_delay (1891, 0), // 8 L2
    make_block_delay (890, 0), // 9 R2

    make_delay (2110), // 10 Loop
    make_quantizer(), // 11
    make_ap (2371, 0.5), // 12 Loop nested allpass 1
    make_ap (1378, 0.2), // 13 Loop nested allpass 2

    make_block_delay (2003, 0), // 14 L3
    make_block_delay (671, 0), // 15 R3

    make_delay (2157 - max_block_size), // 16 Loop
    make_quantizer(), // 17
    make_lp(), // 18
    make_lp(), // 19
    make_ap (2712, 0.5, 22), // 20 Loop nested allpass 1
    make_ap (1783, 0.2), // 21 Loop nested allpass 2

    make_block_delay (max_block_size) // 22 delay block (== blocksz + 1)
  );
}

struct midifex49_spec {
  static constexpr auto values {get_midifex49_spec()};
};
//------------------------------------------------------------------------------
static constexpr auto get_midifex50_spec()
{
  return make_array<stage_data> (
    // outputs first
    make_quantizer(), // 0
    make_ap (147, 0.5), // 1 L diff
    make_ap (43, 0.5), // 2 L diff
    make_ap (55, 0.5), // 3 L diff
    make_ap (249, 0.5), // 4 R diff
    make_ap (48, 0.5), // 5 R diff
    make_ap (21, 0.5), // 6 R diff
    // loop
    make_ap (13, 0.5), // 7
    make_ap (83, 0.5), // 8
    make_ap (116, 0.5), // 9
    make_ap (239, 0.5), // 10
    make_ap (339, 0.5), // 11
    make_ap (481, 0.5), // 12
    make_ap (555, 0.5), // 13
    make_ap (823, 0.5), // 14
    make_ap (999, 0.5), // 15
    make_ap (1100, 0.5), // 16
    make_ap (1347, 0.5), // 17
    make_ap (1563, 0.5), // 18
    make_ap (1841 - 32, 0.5), // 19
    make_block_delay (32), // 20 C1
    make_ap (2001 - 32, 0.5, 67), // 21
    make_block_delay (32), // 22 C2
    make_ap (2083 - 32, 0.5, 127), // 23
    make_lp(), // 24
    make_lp(), // 25
    make_block_delay (max_block_size) // 26 (FB point) (== blocksz + 1)
  );
}

struct midifex50_spec {
  static constexpr auto values {get_midifex50_spec()};
};

//------------------------------------------------------------------------------
static constexpr auto get_dre2000a_spec()
{
  constexpr float g = 0.125f;
  return make_array<stage_data> (
    make_quantizer(), // 0
    make_lp (0.12), // 1
    make_hp (0.99), // 2
    make_parallel_delay (1, 42, -g, 586, g, 1099, -g), // 3 (L)
    make_parallel_delay (1, 105, g, 490, -g, 1290, g), // 4 (R)

    make_comb (5719 - 1, 0.f, 43), // 5
    make_lp(), // 6
    make_lp(), // 7
    make_parallel_delay (
      1, 640, -g, 1494, g, 2199, -g, 3122, g, 4135, -g, 4952, g), // 8

    make_comb (5779 - 1, 0.f, 44), // 9
    make_lp(), // 10
    make_lp(), // 11
    make_parallel_delay (
      1, 902, -g, 1830, g, 2528, -g, 3641, g, 4670, -g, 5432, g), // 12

    make_comb (5905 - 1, 0.f, 42), // 13
    make_lp(), // 14
    make_lp(), // 15
    make_parallel_delay (
      1, 979, -g, 2031, g, 3267, -g, 4104, g, 5103, -g, 5886, g), // 16
    // L
    make_ap (224, -0.7), // 17
    make_delay (311, 311), // 18
    make_ap (227, -0.7), // 19
    make_ap (1343, -0.7), // 20
    // R
    make_ap (224, -0.7), // 21
    make_delay (277, 400), // 22
    make_ap (91, -0.7), // 23
    make_ap (1182, -0.7)); // 24
}

struct dre2000a_spec {
  static constexpr auto values {get_dre2000a_spec()};
};
//------------------------------------------------------------------------------
static constexpr auto get_dre2000b_spec()
{
  constexpr float g = 0.125f;
  return make_array<stage_data> (
    make_quantizer(), // 0
    make_lp (0.15), // 1
    make_hp (0.985), // 2

    make_comb (3821 - 1, 0.f, 54), // 3
    make_lp(), // 4
    make_lp(), // 5
    make_parallel_delay (
      1, 429, -g, 1000, g, 1472, -g, 2088, g, 2765, -g, 3311, g), // 6

    make_comb (4036 - 1, 0.f, 53), // 7
    make_lp(), // 8
    make_lp(), // 9
    make_parallel_delay (
      1, 616, -g, 1225, g, 1691, -g, 2434, g, 3122, -g, 3631, g), // 10

    make_comb (4059 - 1, 0.f, 44.f), // 11
    make_lp(), // 12
    make_lp(), // 13
    make_parallel_delay (
      1, 657, -g, 1359, g, 2184, -g, 2744, g, 3411, -g, 3934, g), // 14

    make_ap (282, -0.7), // 15
    make_ap (343, -0.7), // 16
    // L
    make_delay (311, 311), // 17
    make_ap (233, -0.7), // 18
    make_ap (273, -0.7), // 19
    make_ap (534, -0.7), // 20
    // R
    make_delay (277, 400), // 21
    make_ap (194, -0.7), // 22
    make_ap (426, -0.7), // 23
    make_ap (566, -0.7)); // 24
}

struct dre2000b_spec {
  static constexpr auto values {get_dre2000b_spec()};
};
//------------------------------------------------------------------------------
static constexpr auto get_dre2000c_spec()
{
  constexpr float g = 0.125f;
  return make_array<stage_data> (
    make_quantizer(), // 0
    make_lp (0.18), // 1
    make_hp (0.985), // 2
    make_ap (115, -0.7f), // 3
    make_ap (160, -0.7f), // 4
    make_ap (231, -0.7f), // 5

    make_comb (3794 - 1, 0.f, 23), // 6
    make_lp(), // 7
    make_lp(), // 8
    make_parallel_delay (2, 20 + 13, -g, 959, g, 1817, -g, 2855, g), // 9 R/L

    make_comb (3838 - 1, 0.f, 24), // 10
    make_lp(), // 11
    make_lp(), // 12
    make_parallel_delay (2, 339, -g, 1309, g, 2271, -g, 3221, g), // 13 R/L

    make_comb (3861 - 1, 0.f, 25), // 14
    make_lp(), // 15
    make_lp(), // 16
    make_parallel_delay (2, 130, -g, 1104, g, 2065, -g, 3391, g), // 17 R/L

    make_comb (3894 - 1, 0.f, 26), // 18
    make_lp(), // 19
    make_lp(), // 20
    make_parallel_delay (2, 499, -g, 1445, g, 2071, -g, 2885, g), // 21 R/L
    // L
    make_ap (140, -0.7), // 22
    make_ap (160, -0.7), // 23
    make_delay (212, 237), // 24
    // R
    make_ap (160, -0.7), // 25
    make_ap (227, -0.7), // 26
    make_delay (213, 126) // 27
  );
}

struct dre2000c_spec {
  static constexpr auto values {get_dre2000c_spec()};
};
//------------------------------------------------------------------------------
static constexpr auto get_dre2000d_spec()
{
  constexpr float g = 0.125f;
  return make_array<stage_data> (
    make_quantizer(), // 0
    make_lp (0.185), // 1
    make_hp (0.99), // 2
    make_ap (115, -0.7f), // 3
    make_ap (214, -0.7f), // 4
    make_ap (365, -0.7f), // 5

    make_comb (5652 - 1, 0.f, 57), // 6
    make_lp(), // 7
    make_lp(), // 8
    make_parallel_delay (2, 745, -g, 2164, g, 2922, -g, 5005, g), // 9 R/L

    make_comb (5679 - 1, 0.f), // 10
    make_lp(), // 11
    make_lp(), // 12
    make_parallel_delay (2, 625, -g, 2080, g, 3523, -g, 4948, g), // 13 R/L

    make_comb (5689 - 1, 0.f), // 14
    make_lp(), // 15
    make_lp(), // 16
    make_parallel_delay (2, 11 + 21, -g, 1435, g, 2769, -g, 4279, g), // 17 R/L

    make_comb (5711 - 1, 0.f, 57), // 18
    make_lp(), // 19
    make_lp(), // 20
    make_parallel_delay (2, 59, -g, 1557, g, 3539, -g, 4304, g), // 21 R/L
    // L
    make_ap (214, -0.7, 30), // 22
    make_ap (469, -0.7), // 23
    make_delay (234, 337), // 24
    // R
    make_ap (245, -0.7, 29), // 25
    make_ap (426, -0.7), // 26
    make_delay (434, 247) // 27
  );
}

struct dre2000d_spec {
  static constexpr auto values {get_dre2000d_spec()};
};

//------------------------------------------------------------------------------
static constexpr auto get_rev5_l_hall_spec()
{
  constexpr float g = 0.125f;
  return make_array<stage_data> (
    make_quantizer(), // 0
    make_lp (0.32), // 1
    make_hp (0.99), // 2

    make_quantizer(), // 3
    make_ap (1257, 0.75f), // 4
    make_ap (978, 0.75f), // 5

    // rev main
    make_comb (2076, 0., 43), // 6
    make_lp(), // 7
    make_lp(), // 8
    make_parallel_delay (2, 949, g, 837, g), // 9
    make_comb (2894, 0., 43), // 10
    make_lp(), // 11
    make_lp(), // 12
    make_parallel_delay (2, 1330, g, 1412, g), // 13
    make_comb (3295, 0.), // 14
    make_lp(), // 15
    make_lp(), // 16
    make_parallel_delay (2, 2393, g, 2704, g), // 17
    make_comb (3919, 0.), // 18
    make_lp(), // 19
    make_lp(), // 20
    make_parallel_delay (2, 3263, g, 3011, g), // 21
    make_comb (4570, 0.), // 22
    make_lp(), // 23
    make_lp(), // 24
    make_parallel_delay (2, 3667, g, 3667, g), // 25
    // sub reverb cond
    make_variable_delay (533, 533 + 2411), // 26
    make_ap (1212, 0.75f), // 27
    make_ap (960, 0.75f), // 28
    // sub reverb
    make_comb (3169, 0., 44), // 29
    make_lp(), // 30
    make_lp(), // 31
    make_parallel_delay (2, 1757, g, 1644, g), // 32
    make_comb (3753, 0., 44), // 33
    make_lp(), // 34
    make_lp(), // 35
    make_parallel_delay (2, 1900, g, 1981, g), // 36
    make_comb (4280, 0.), // 37
    make_lp(), // 38
    make_lp(), // 39
    make_parallel_delay (2, 2838, g, 3148, g), // 40
    make_comb (4491, 0.), // 41
    make_lp(), // 42
    make_lp(), // 43
    make_parallel_delay (2, 3798, g, 3545, g), // 44
    make_comb (5091, 0.), // 45
    make_lp(), // 46
    make_lp(), // 47
    make_parallel_delay (2, 4298, g, 4298, g), // 48
    make_ap (682, 0.75f), // 49
    make_ap (830, 0.75f), // 50
    make_ap (695, 0.75f), // 51
    make_ap (844, 0.75f), // 52
    make_parallel_delay (
      2,
      18 + 14,
      0.95f,
      18 + 14,
      0.95f,
      895,
      0.95f,
      695,
      0.95f,
      1238,
      0.75f,
      1487,
      0.75f,
      1870,
      0.75f,
      1936,
      0.75f,
      2724,
      0.5f,
      2635,
      0.5f,
      3365,
      0.25f,
      3609,
      0.25f) // 53 L/R
  );
}

struct rev5_l_hall_spec {
  static constexpr auto values {get_rev5_l_hall_spec()};
};

}} // namespace detail::lofiverb
//------------------------------------------------------------------------------
// A reverb using 16-bit fixed-point arithmetic on the main loop. One design
// criteria has been for it to be extremely CPU friendly.
// Notice: __fp16 (half-precision floating point storage) could achieve the
// same memory savings with more dynamic range. This is still kept as a 16-bit
// fixed-point reverb.
class lofiverb {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct algorithm_tag {};
  void set (algorithm_tag, int v)
  {
    u32 curr = _param.mode / 2;
    if (v == curr) {
      return;
    }
    _param.mode = (_param.mode % 2) | (v * 2);
    update_operation_mode();
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (algorithm_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Artv Abyss",
        "Artv Small",
        "Acreil Midifex 49",
        "Acreil Midifex 50",
        "Acreil Dre-2000 A",
        "Acreil Dre-2000 B",
        "Acreil Dre-2000 C",
        "Acreil Dre-2000 D",
        "Acreil REV5 L Hall"
#ifdef LOFIVERB_ADD_DEBUG_ALGO
        ,
        "Debug "
#endif
        ),
      64);
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void set (mode_tag, int v)
  {
    u32 curr = _param.mode & 1;
    if (v == curr) {
      return;
    }
    _param.mode = (_param.mode & ~(1u)) | (v & 1);
    update_operation_mode();
  }

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      1, make_cstr_array ("16-bit float point", "16-bit fixed point"), 16);
  }
  //----------------------------------------------------------------------------
  struct clock_tag {};
  void set (clock_tag, int v)
  {
    if (v == _param.srateid) {
      return;
    }
    _param.srateid = v;
    update_operation_mode();
  }

  static constexpr auto get_parameter (clock_tag)
  {
    return choice_param (
      3,
      make_cstr_array (
        "Downclocked 3",
        "Downclocked 2",
        "Downclocked 1",
        "Base Clock",
        "Overclocked 1",
        "Overclocked 2",
        "Overclocked 3"));
  }
  //----------------------------------------------------------------------------
  struct decay_tag {};
  void set (decay_tag, float v)
  {
    _param_smooth.target().decay = v * 0.01f * fixpt_max_flt;
  }

  static constexpr auto get_parameter (decay_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct mod_tag {};
  void set (mod_tag, float v)
  {
    v *= 0.01f * fixpt_max_flt;
    if (v == _param_smooth.target().mod) {
      return;
    }
    _param_smooth.target().mod = v;
    update_mod();
  }

  static constexpr auto get_parameter (mod_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct character_tag {};
  void set (character_tag, float v)
  {
    _param_smooth.target().character = v * 0.01f * fixpt_max_flt;
  }

  static constexpr auto get_parameter (character_tag)
  {
    return float_param ("%", 0., 100., 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lf_amt_tag {};
  void set (lf_amt_tag, float v) { _param.lf_amt = v * 0.01f * fixpt_max_flt; }

  static constexpr auto get_parameter (lf_amt_tag)
  {
    return float_param ("%", 0., 100., 50., 0.01);
  }
  //----------------------------------------------------------------------------
  struct hf_amt_tag {};
  void set (hf_amt_tag, float v) { _param.hf_amt = v * 0.01f * fixpt_max_flt; }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (hf_amt_tag)
  {
    return float_param ("%", 0., 100., 50., 0.01);
  }
  //----------------------------------------------------------------------------
  static constexpr uint max_predelay_qb = 4;
  struct predelay_tag {};
  void set (predelay_tag, float v) { _param.predelay = v; }

  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("quarters", 0., max_predelay_qb, 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct clip_level_tag {};
  void set (clip_level_tag, float v)
  {
    _param.gain = db_to_gain (v) * clip_value();
  }

  static constexpr auto get_parameter (clip_level_tag)
  {
    return float_param ("dBFS", -36., 36., 0., 0.1, 0.5, true);
  }
  //----------------------------------------------------------------------------
  struct stereo_tag {};
  void set (stereo_tag, float v) { _param_smooth.target().stereo = v * 0.01f; }

  static constexpr auto get_parameter (stereo_tag)
  {
    return float_param ("%", -100., 100., 100., 0.01);
  }
  //----------------------------------------------------------------------------
  struct ducking_threshold_tag {};
  void set (ducking_threshold_tag, float v)
  {
    if (v == _param.ducking_threshold) {
      return;
    }
    _param.ducking_threshold = v;
    _ducker.set_threshold (vec_set<2> (v));
  }

  static constexpr auto get_parameter (ducking_threshold_tag)
  {
    return float_param ("dB", -60.f, 0.0f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct ducking_speed_tag {};
  void set (ducking_speed_tag, float v)
  {
    if (v == _param.ducking_speed) {
      return;
    }
    _param.ducking_speed = v;
    update_ducker();
  }

  static constexpr auto get_parameter (ducking_speed_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _pc                  = &pc;
    _n_processed_samples = 0;
    _srate               = 0; // ensure triggering a resample reset
    _param.mode
      = (decltype (_param.mode)) -1ull; // trigger a mode and resampler reset
    _param.srateid = 3;

    _param_smooth.reset (_t_spl);
    // set defaults
    mp11::mp_for_each<parameters> ([&] (auto param) {
      set (param, get_parameter (param).min);
      if constexpr (!is_choice<decltype (get_parameter (param))>) {
        set (param, get_parameter (param).max);
      }
      else {
        // max might be not yet impl.
        set (param, get_parameter (param).min + 1);
      }
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
#if MOD_DBG_DENORMALS
    feenableexcept (FE_INVALID);
#endif
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    _resampler.process (outs, ins, samples, [=] (auto io) {
      if ((_param.mode % 2) == 0) {
        process_float (io);
      }
      else {
        process_fixed_point (io);
      }
    });
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    algorithm_tag, // keep first! important for reset to work
    mode_tag,
    clock_tag,
    character_tag,
    lf_amt_tag,
    decay_tag,
    predelay_tag,
    hf_amt_tag,
    clip_level_tag,
    mod_tag,
    stereo_tag,
    ducking_threshold_tag,
    ducking_speed_tag>;
  //----------------------------------------------------------------------------
private:
  // fixpt_t truncates when dropping fractional bits (leaky).
  // fixpt_tr rounds to the nearest (never reaches full zero).
  //
  // Using fixpt_t as default, with fixpt_tr at some points compensate the
  // truncating leakage.
  using fixpt_t  = detail::lofiverb::fixpt_t;
  using fixpt_tr = detail::lofiverb::fixpt_tr;
  using float16  = detail::lofiverb::float16;

  static constexpr auto blank     = detail::lofiverb::defaulted;
  static constexpr auto add_to    = detail::lofiverb::add_to;
  static constexpr auto overwrite = detail::lofiverb::overwrite;
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = detail::lofiverb::max_block_size;
  static constexpr uint n_channels     = 2;
  // this is using 16 bits fixed-point arithmetic, positive values can't
  // represent one, so instead of correcting everywhere the parameters are
  // scaled instead to never reach 1.
  static constexpr auto fixpt_max_flt = fixpt_t::max_float();
  //----------------------------------------------------------------------------
  static float clip_value() { return std::min (fixpt_max_flt, float16::max()); }
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr auto one()
  {
    if constexpr (std::is_floating_point_v<T>) {
      return fixpt_max_flt;
    }
    else {
      return fixpt_t::max();
    }
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters;
  struct smoothed_parameters;
  //----------------------------------------------------------------------------
  void process_fixed_point (xspan<std::array<float, 2>> io)
  {
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    array2d<fixpt_t, 2, max_block_size> wet;
    std::array<f32_x2, max_block_size>  ducker_gain;
    loop_parameters                     pars;

    auto pre_gain = (1.f / _param.gain);

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto wetv      = vec_from_array (io[i]) * pre_gain;
      wetv           = vec_clamp (wetv, -clip_value(), clip_value());
      ducker_gain[i] = _ducker.tick (wetv);
      io[i]          = vec_to_array (wetv);
      _param_smooth.tick();
      pars.stereo[i] = _param_smooth.get().stereo;
      pars.decay[i].load_float (_param_smooth.get().decay);
      pars.character[i].load_float (_param_smooth.get().character);
      pars.mod[i].load_float (_param_smooth.get().mod);
    }
    // predelay + int conversion
    if (_param.predelay != 0) {
      uint predelay_spls = _1_4beat_spls * _param.predelay;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        wet[i][0].load_float (_predelay.get (predelay_spls, 0));
        wet[i][1].load_float (_predelay.get (predelay_spls, 1));
        _predelay.push (io[i]);
      }
    }
    else {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        wet[i][0].load_float (io[i][0]);
        wet[i][1].load_float (io[i][1]);
      }
    }
    // main loop
    switch (_param.mode) {
#ifdef LOFIVERB_ADD_DEBUG_ALGO
    case mode::debug_algo:
      process_debug_algo (xspan {wet.data(), io.size()}, pars);
      break;
#endif
    case mode::abyss:
      process_abyss (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::small_space:
      process_small_space (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::midifex49:
      process_midifex49 (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::midifex50:
      process_midifex50 (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::dre2000a:
      process_dre2000a (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::dre2000b:
      process_dre2000b (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::dre2000c:
      process_dre2000c (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::dre2000d:
      process_dre2000d (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::rev5_l_hall:
      process_rev5_l_hall (xspan {wet.data(), io.size()}, pars);
      break;
    default:
      assert (false);
    }
    // float conversion
    auto postgain = _param.gain;
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = wet[i][0].to_floatp();
      auto r = wet[i][1].to_floatp();
      l *= ducker_gain[i][0] * postgain;
      r *= ducker_gain[i][1] * postgain;
      l = r * (1 - abs (pars.stereo[i])) + l * abs (pars.stereo[i]);
      if (pars.stereo[i] < 0) {
        io[i][0] = r;
        io[i][1] = l;
      }
      else {
        io[i][0] = l;
        io[i][1] = r;
      }
    }
  }
  //----------------------------------------------------------------------------
  void process_float (xspan<std::array<float, 2>> io)
  {
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    std::array<f32_x2, max_block_size> ducker_gain;
    loop_parameters_flt                pars;

    auto pre_gain = 1.f / _param.gain;

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto wet       = vec_from_array (io[i]) * pre_gain;
      wet            = vec_clamp (wet, -clip_value(), clip_value());
      ducker_gain[i] = _ducker.tick (wet);
      io[i]          = vec_to_array (wet);
      _param_smooth.tick();
      pars.stereo[i]    = _param_smooth.get().stereo;
      pars.decay[i]     = _param_smooth.get().decay;
      pars.character[i] = _param_smooth.get().character;
      pars.mod[i]       = _param_smooth.get().mod;
    }
    // predelay
    if (_param.predelay != 0) {
      uint predelay_spls = _1_4beat_spls * _param.predelay;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        // TODO: block fetch?
        std::array<float, 2> spl;
        spl[0] = _predelay.get (predelay_spls, 0);
        spl[1] = _predelay.get (predelay_spls, 1);
        _predelay.push (xspan {io[i]});
        io[i] = spl;
      }
    }
    // main loop
    switch (_param.mode) {
#ifdef LOFIVERB_ADD_DEBUG_ALGO
    case mode::debug_algo_flt:
      process_debug_algo (io, pars);
      break;
#endif
    case mode::abyss_flt:
      process_abyss (io, pars);
      break;
    case mode::small_space_flt:
      process_small_space (io, pars);
      break;
    case mode::midifex49_flt:
      process_midifex49 (io, pars);
      break;
    case mode::midifex50_flt:
      process_midifex50 (io, pars);
      break;
    case mode::dre2000a_flt:
      process_dre2000a (io, pars);
      break;
    case mode::dre2000b_flt:
      process_dre2000b (io, pars);
      break;
    case mode::dre2000c_flt:
      process_dre2000c (io, pars);
      break;
    case mode::dre2000d_flt:
      process_dre2000d (io, pars);
      break;
    case mode::rev5_l_hall_flt:
      process_rev5_l_hall (io, pars);
      break;

    default:
      assert (false);
    }
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = io[i][0] * ducker_gain[i][0] * _param.gain;
      auto r = io[i][1] * ducker_gain[i][1] * _param.gain;
      l      = r * (1 - abs (pars.stereo[i])) + l * abs (pars.stereo[i]);
      if (pars.stereo[i] < 0) {
        io[i][0] = r;
        io[i][1] = l;
      }
      else {
        io[i][0] = l;
        io[i][1] = r;
      }
    }
  }
//----------------------------------------------------------------------------
#ifdef LOFIVERB_ADD_DEBUG_ALGO
  template <class T, class Params>
  void process_debug_algo (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<debug_algo_type> (_modes);

    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr    in_arr;
    arr_fb loop_arr;

    auto in   = xspan {in_arr.data(), io.size()};
    auto loop = xspan {loop_arr.data(), io.size() + 1};

    span_visit (in, [&] (auto& v, uint i) {
      v = (T) ((io[i][0] + io[i][1]) * 0.25_r);
    });
    rev.fetch_block<3> (loop, 1);
    span_visit (io, [&] (auto& v, uint i) { v[0] = v[1] = loop[i + 1]; });
    loop.cut_tail (1);
    span_visit (loop, [&] (auto& v, uint i) {
      v = (T) (in[i] + v * (0.7_r + par.decay[i] * 0.28_r));
    });
    rev.run<0, 1, 2> (loop);
    rev.push<3> (loop.to_const()); // feedback point
  }
#endif
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_abyss (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<abyss_type> (_modes);

    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr late_in_arr;
    arr lfo1;
    arr lfo2;
    arr lfo3;
    arr lfo4;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = (T) ((io[i][0] + io[i][1]) * 0.5_r);
      auto mod   = (T) (0.25_r + (1_r - par.mod[i]) * 0.75_r);
      // ER + late lfo
      auto lfo = tick_lfo<T>();
      lfo1[i]  = T {lfo[0]};
      lfo2[i]  = (T) (T {lfo[1]} * (0.5_r + par.character[i] * 0.5_r));
      lfo3[i]  = (T) (T {lfo[2]} * mod);
      lfo4[i]  = (T) (T {lfo[3]} * mod);

      // decay fixup
      auto decay   = (T) (one<T>() - par.decay[i]);
      decay        = (T) (one<T>() - decay * decay);
      par.decay[i] = (T) (0.6_r + decay * 0.38_r);
    }
    // damp -----------------------------------
    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.55f + _param.lf_amt * 0.4f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.45f + _param.hf_amt * 0.4f);

    // diffusion -----------------------------
    rev.run<0> (late_in);
    rev.run<1> (late_in);
    rev.run<2> (late_in);
    rev.run<3> (late_in);

    // ER (first reverb, not exactly ER at the end...) -----------------------
    arr    early1_arr;
    arr    early1b_arr;
    arr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    rev.fetch_block<9> (er2, 1); // feedback, fetching block + 1 samples

    span_visit (er1, [&] (auto& v, uint i) {
      v = (T) (late_in[i] * 0.5_r + er2[i]);
    });
    er2.cut_head (1); // drop feedback sample from previous block

    rev.run<4> (er1, xspan {lfo2});
    rev.run<5, 6> (er1, flo, glo, fhi, ghi);
    xspan_memcpy (er1b, er1);
    rev.run<7> (er1b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<8> (er1b, xspan {lfo1});
    rev.push<9> (er1b.to_const()); // feedback point

    // Late -----------------------------
    arr    late_arr;
    arr_fb l_arr;
    arr_fb r_arr;
    arr    g_arr;

    auto late = xspan {late_arr.data(), io.size()};
    auto g    = xspan {g_arr.data(), io.size()};
    auto l    = xspan {l_arr.data(), io.size() + 1}; // +1: Feedback on head
    auto r    = xspan {r_arr.data(), io.size() + 1}; // +1: Feedback on head

    // feedback handling
    rev.fetch_block<19> (l, 1); // feedback, fetching block + 1 samples
    rev.fetch_block<30> (r, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      auto loopsig = late_in[i] * 0.5_r + r[i];
      auto er_sig  = (er1[i] + er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (one<T>() - er_amt) + er_sig * er_amt);
      g[i] = (T) (0.618_r + par.character[i] * ((0.707_r - 0.618_r) * 2_r));
    }
    r.cut_head (1); // drop feedback sample from previous block

    rev.run<10> (late, xspan {lfo3}, g);
    rev.run<11> (late);
    rev.run<12, 13> (late, flo, glo, fhi, ghi);
    rev.run<14> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<15> (late, blank, [g] (uint i) { return -g[i]; });
    rev.run<16> (late);
    rev.run<17> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<18> (late);
    rev.push<19> (late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      auto loopsig = late_in[i] * 0.5_r + l[i];
      auto er_sig  = (er1[i] - er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (one<T>() - er_amt) + er_sig * er_amt);
    }
    l.cut_head (1); // drop feedback sample from previous block
    rev.run<20> (late, xspan {lfo4}, g);
    rev.run<21> (late);
    rev.run<22, 23> (late, flo, glo, fhi, ghi);
    rev.run<24> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<25> (late, blank, [g] (uint i) { return -g[i]; });
    rev.run<26> (late);
    rev.run<27> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<28> (late);
    rev.run<29> (late);
    rev.push<30> (late.to_const()); // feedback point

    // Mixdown
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto er_amt     = par.character[i] * 0.2_r;
      auto e_l        = (-er1[i] * 0.66_r - er2[i] * 0.34_r) * er_amt;
      auto e_r        = (-er1[i] * 0.66_r + er2[i] * 0.34_r) * er_amt;
      auto direct_amt = one<T>() - er_amt;
      io[i][0]        = (T) (l[i] * direct_amt + e_l);
      io[i][1]        = (T) (r[i] * direct_amt + e_r);
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_small_space (xspan<std::array<T, 2>> io, Params& par)
  {
    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    auto& rev = std::get<small_space_type> (_modes);

    arr_fb a_mem, d_mem;
    arr    b_mem, c_mem, tmp1, tmp2, tmp3;
    arr    lfo1, lfo2;

    xspan a {a_mem.data(), io.size() + 1};
    xspan b {b_mem.data(), io.size()};
    xspan c {c_mem.data(), io.size()};
    xspan d {d_mem.data(), io.size() + 1};

    // fetch a block from the output (a and d)
    rev.fetch_block<47> (a, 1); // feedback, fetching block + 1 samples
    rev.fetch_block<67> (d, 1); // feedback, fetching block + 1 samples

    xspan_memcpy (b, a.advanced (1)); // L out on b
    a.cut_tail (1); // fb samples on a.

    xspan_memcpy (c, d.advanced (1)); // R out on c
    d.cut_tail (1); // fb samples on d.

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // invert character first
      par.character[i] = (T) (one<T>() - par.character[i]);
      tmp1[i]          = (T) (0.4_r + 0.3_r * par.character[i]);
      tmp2[i]          = (T) (-0.6_r * par.character[i]);
      tmp3[i]          = (T) (0.6_r * par.character[i]);
      // decay fixup
      auto d       = one<T>() - par.decay[i];
      par.decay[i] = (T) (0.9_r - d * d * 0.2_r);
    }

    // output preconditioning
    rev.run<0> (b, 0, tmp1);
    rev.run<1> (b, 0, tmp2);
    rev.run<2> (b, 0, tmp2);
    rev.run<3> (b, 0, tmp3);

    rev.run<4> (c, 0, tmp1);
    rev.run<5> (c, 0, tmp2);
    rev.run<6> (c, 0, tmp2);
    rev.run<7> (c, 0, tmp3);

    xspan i1 {tmp1.data(), io.size()};
    xspan i2 {tmp2.data(), io.size()};
    // output write + preparations
    span_visit (io, [&] (auto& spls, uint i) {
      i1[i]   = (T) ((spls[0] + spls[1]) * 0.25_r); // m
      i2[i]   = (T) ((spls[0] - spls[1]) * 0.25_r); // s
      spls[0] = b[i];
      spls[1] = c[i];
      b[i]    = (T) (0.8_r + 0.1_r * par.character[i]); // character 1 on b
      c[i]    = (T) (0.8_r * par.character[i]); // character 2 on c
    });

    // input preconditioning
    rev.run<8> (i1, 0, b);
    rev.run<9> (i1, 0, c);
    rev.run<10> (i1, 0, c);

    rev.run<11> (i2, 0, b);
    rev.run<12> (i2, 0, c);
    rev.run<13> (i2, 0, c);

    // fetch b and d feedback (a and d are already ready) 1 block only
    rev.fetch_block<54> (b, 1); // feedback
    rev.fetch_block<61> (c, 1); // feedback

    // quantized decay
    rev.run<14> (a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<15> (b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<16> (c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<17> (d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // sum inputs to feedback
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // a single hadamard iteration can duplicate the range on one
      // of the channels, so halving once more.
      b[i] = (T) (b[i] + tmp1[i] * 0.5_r);
      d[i] = (T) (d[i] + tmp2[i] * 0.5_r);
    }
    rev.hadamard4 (
      make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto lfo = tick_lfo<T>();
      lfo1[i]  = T {lfo[0]};
      lfo2[i]  = T {lfo[1]};
      // lfo3 = -lfo1, lfo4 = -lfo2
      // mod fixup
      par.mod[i] = (T) (one<T>() - par.mod[i]);
      par.mod[i] = (T) (one<T>() - (par.mod[i] * par.mod[i]));
      // preparing modlations for block a
      tmp1[i] = (T) (0.1_r + 0.1_r * lfo1[i]);
    }
    // channel a block 1
    rev.run<18, 19, 20> (a, blank, tmp1);
    rev.run<21> (a);

    // channel b block 1
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (-0.1_r - 0.1_r * -lfo1[i]); // lfo3 = -lfo1
    }
    rev.run<22, 23, 24> (b, 0, tmp1);
    rev.run<25> (b);

    // channel c block 1
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (0.1_r + 0.1_r * -lfo2[i]); // lfo4 = -lfo2
    }
    rev.run<26, 27, 28> (c, 0, tmp1);
    rev.run<29> (c);

    // channel d block 1
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (0.1_r - 0.1_r * lfo2[i]);
    }
    rev.run<30, 31, 32> (d, 0, tmp1);
    rev.run<33> (d);

    // quantized decay
    rev.run<34> (a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<35> (b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<36> (c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<37> (d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // block 2
    rev.hadamard4 (
      make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // swaps.
    a = xspan {d_mem.data(), io.size()};
    b = xspan {c_mem.data(), io.size()};
    c = xspan {a_mem.data(), io.size()};
    d = xspan {b_mem.data(), io.size()};

    // damp -----------------------------------
    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.87f + _param.lf_amt * 0.13f);
    T fhi = load_float<T> (0.75f - _param.hf_amt * _param.hf_amt * 0.55f);
    T ghi = load_float<T> (0.7f + _param.hf_amt * 0.24f);

    // channel a block 2
    rev.run<38, 39> (a, flo, glo, fhi, ghi);
    rev.run<40, 41, 42> (a);
    rev.run<43> (a);
    rev.run<44> (a);
    rev.run<45> (a);
    rev.run<46> (a);
    rev.push<47> (a.to_const());

    // channel b block 2
    rev.run<48, 49> (b, flo, glo, fhi, ghi);
    rev.run<50, 51, 52> (b);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (lfo1[i] * par.mod[i]);
      tmp2[i] = (T) (0.4_r + 0.4_r * lfo2[i] * par.mod[i]);
    }
    rev.run<53> (b, tmp1, tmp2);
    rev.push<54> (b.to_const());

    // channel c block 2
    rev.run<55, 56> (c, flo, glo, fhi, ghi);
    rev.run<57, 58, 59> (c);
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (-tmp1[i]); // lfo3 = -lfo1.
      tmp2[i] = (T) (0.6_r - 0.4_r * -lfo2[i] * par.mod[i]); // lfo4 = -lfo2
    }
    rev.run<60> (c, tmp1, tmp2);
    rev.push<61> (c.to_const());

    // channel d block 2
    rev.run<62, 63> (d, flo, glo, fhi, ghi);
    rev.run<64, 65, 66> (d);
    rev.push<67> (d.to_const());
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_midifex49 (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<midifex49_type> (_modes);

    using arr = block_arr<T>;

    arr   tmp_arr;
    xspan tmp {tmp_arr.data(), io.size()};
    arr   sig_arr;
    xspan sig {sig_arr.data(), io.size()};
    arr   l_arr;
    xspan l {l_arr.data(), io.size()};
    arr   r_arr;
    xspan r {r_arr.data(), io.size()};
    arr   lfo1_arr;
    xspan lfo1 {lfo1_arr.data(), io.size()};
    arr   lfo2_arr;
    xspan lfo2 {lfo2_arr.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      // Midside signal
      sig[i] = (T) (((spl[0] + spl[1]) * 0.25_r)); // gain = 0.5 + feedback = 1
      // LFO
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      // decay fixup
      auto decay   = (T) (one<T>() - par.decay[i]);
      decay        = (T) (one<T>() - decay * decay);
      par.decay[i] = (T) (0.1_r + decay * 0.8375_r);
    });

    rev.run<0> (sig);
    rev.run<1> (sig);
    rev.run<2> (sig);
    rev.run<3> (sig);

    // feedback handling, fetching the block with a negative offset of 1
    rev.fetch_block<22> (tmp, 1);
    span_visit (sig, [&] (auto& v, uint i) { v = (T) (v + tmp[i]); });

    // 1st output point for L and R signal
    xspan_memcpy (l, sig); // delay LT a block -> might overlap, requires copy
    rev.run<4> (l);
    rev.fetch_block<5> (r, 0); // delay GT a block will never overlap
    rev.push<5> (sig.to_const());

    // continuing the loop
    rev.run<6> (sig);
    rev.run<7> (sig, lfo1, blank);

    // 2nd output point for L and R signal
    rev.fetch_block<8> (tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push<8> (sig.to_const());
    rev.fetch_block<9> (tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push<9> (sig.to_const());

    // continuing the loop
    rev.run<10> (sig);
    rev.run<11> (sig, [&] (auto v, uint i) { return v * par.decay[i]; });
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.14_r);
    });
    rev.run<12, 13> (sig, blank, blank, blank, tmp);

    // 3rd output point for L and R signal
    rev.fetch_block<14> (tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push<14> (sig.to_const());
    rev.fetch_block<15> (tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push<15> (sig.to_const());

    // continuing the loop
    rev.run<16> (sig);
    rev.run<17> (sig, [&] (auto v, uint i) { return v * par.decay[i]; });

    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.35f + _param.lf_amt * 0.6f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.35f + _param.hf_amt * 0.5f);

    rev.run<18, 19> (sig, flo, glo, fhi, ghi);
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.2_r);
    });
    rev.run<20, 21> (sig, lfo2, blank, blank, tmp);
    // push to delay feedback
    rev.push<22> (sig.to_const());

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_midifex50 (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<midifex50_type> (_modes);

    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr_fb tmp1;
    xspan  loop {tmp1.data(), io.size() + 1};
    arr    l_arr;
    xspan  l {l_arr.data(), io.size()};
    arr    r_arr;
    xspan  r {r_arr.data(), io.size()};

    // xspan_memcpy (loop, io);
    rev.fetch_block<26> (loop, 1);
    xspan_memcpy (r, loop.advanced (1));
    xspan_memcpy (l, loop.advanced (1));

    loop.cut_tail (1); // feedback samples
    rev.run<0> (loop, [&] (T fb_spl, uint i) {
      // decay fixup
      auto decay = (T) (one<T>() - par.decay[i]);
      decay      = (T) (one<T>() - decay * decay);
      decay      = (T) - (0.3_r + decay * 0.45_r);
      return (fb_spl * decay) + ((io[i][0] + io[i][1]) * 0.25_r); // gain = 1
    }); // feedback + input summing with quantizer

    rev.run<1> (l);
    rev.run<2> (l);
    rev.run<3> (l);

    rev.run<4> (r);
    rev.run<5> (r);
    rev.run<6> (r);

    rev.run<7> (loop);
    rev.run<8> (loop);
    rev.run<9> (loop);
    rev.run<10> (loop);
    rev.run<11> (loop);
    rev.run<12> (loop);
    rev.run<13> (loop);
    rev.run<14> (loop);
    rev.run<15> (loop);
    rev.run<16> (loop);
    rev.run<17> (loop);
    rev.run<18> (loop);
    rev.run<19> (loop);

    arr tmp2;
    arr tmp3;

    xspan c1 {tmp2.data(), loop.size()};
    xspan lfo2 {tmp3.data(), loop.size()};
    rev.fetch_block<20> (c1, 0);
    rev.push<20> (loop.to_const());
    auto lfo1 = loop; // data inserted already (lfo1 -> tmp1)
    loop      = c1; // avoid a copy. loop -> tmp2
    span_visit (l, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = one<T>() - c;
      v      = (T) (k * v + loop[i] * c); // L done
      // unrelated but done here to skip one iteration
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
    });
    rev.run<21> (loop, xspan {lfo1}, blank); // tmp1 free
    xspan c2 {tmp1.data(), loop.size()};
    rev.fetch_block<22> (c2, 0);
    rev.push<22> (loop.to_const());
    loop = c2; // avoid copy. tmp2 (loop) free.
    span_visit (r, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = one<T>() - c;
      v      = (T) (k * v + loop[i] * c); // R done
    });
    // outputs done, dump now that they have been recently touched
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
    rev.run<23> (loop, xspan {lfo2}, blank); // tmp3 free

    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.35f + _param.lf_amt * 0.6f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.35f + _param.hf_amt * 0.5f);

    rev.run<24, 25> (loop, flo, glo, fhi, ghi);
    rev.push<26> (loop.to_const());
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_dre2000a (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<dre2000a_type> (_modes);

    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, tmp1, tmp2, tank;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (l, in.to_const());
    xspan_memdump (r.data(), l);
    rev.run<3> (l, overwrite);
    rev.run<4> (r, overwrite);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 5, 9, 13> (0.9f + dec2 * 19.1f, _srate);
    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.8f + _param.lf_amt * 0.18f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.4f + dec2 * 0.15f + _param.hf_amt * 0.35f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<5> (comb_fb, lfo1, -gains[0]);
    rev.run<6, 7> (comb_fb, flo, glo, fhi, ghi);
    rev.push<5> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<8> (comb_fb.to_const(), overwrite, tank);

    rev.fetch_block<9> (comb_fb, lfo2, -gains[1]);
    rev.run<10, 11> (comb_fb, flo, glo, fhi, ghi);
    rev.push<9> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<12> (comb_fb.to_const(), add_to, tank);

    rev.fetch_block<13> (comb_fb, lfo3, -gains[2]);
    rev.run<14, 15> (comb_fb, flo, glo, fhi, ghi);
    rev.push<13> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<16> (comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      l[i]     = (T) (l[i] + tank[i]);
      r[i]     = (T) (r[i] + tank[i]);
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    auto stank = xspan {tank.data(), io.size()}.to_const();
    rev.run<17> (l, stank);
    xspan er {tmp2.data(), io.size()};
    rev.run<18> (er, in.to_const(), par.character);
    crossfade (l, er, eramt);
    rev.run<19> (l);
    rev.run<20> (l);

    rev.run<21> (r, stank);
    rev.run<22> (er, in.to_const(), par.character);
    crossfade (r, er, eramt);
    rev.run<23> (r);
    rev.run<24> (r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_dre2000b (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<dre2000b_type> (_modes);

    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, tmp1, tmp2, tank;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 3, 7, 11> (0.25f + dec2 * 10.f, _srate);
    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.88f + _param.lf_amt * 0.1f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.25f + dec2 * 0.4f + _param.hf_amt * 0.3f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<3> (comb_fb, lfo1, -gains[0]);
    rev.run<4, 5> (comb_fb, flo, glo, fhi, ghi);
    rev.push<3> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<6> (comb_fb.to_const(), overwrite, tank);

    rev.fetch_block<7> (comb_fb, lfo2, -gains[1]);
    rev.run<8, 9> (comb_fb, flo, glo, fhi, ghi);
    rev.push<7> (comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (T) (s + in[i]); });
    rev.run<10> (comb_fb.to_const(), add_to, tank);

    rev.fetch_block<11> (comb_fb, lfo3, -gains[2]);
    rev.run<12, 13> (comb_fb, flo, glo, fhi, ghi);
    rev.push<11> (comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (T) (s + in[i]); });
    rev.run<14> (comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    xspan stank {tank.data(), io.size()};
    rev.run<15> (stank);
    rev.run<16> (l, stank.to_const());
    xspan_memdump (r.data(), l);

    xspan er {tmp2.data(), io.size()};
    rev.run<17> (er, in.to_const(), par.character);
    crossfade (l, er, eramt);
    rev.run<18> (l);
    rev.run<19> (l);
    rev.run<20> (l);

    rev.run<21> (er, in.to_const(), par.character);
    crossfade (r, er, eramt);
    rev.run<22> (r);
    rev.run<23> (r);
    rev.run<24> (r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_dre2000c (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<dre2000c_type> (_modes);

    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, lfo4, tmp1;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);
    rev.run<3> (in);
    rev.run<4> (in);
    rev.run<5> (in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 6, 10, 14, 18> (0.25f + dec2 * 20.f, _srate);
    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.88f + _param.lf_amt * 0.1f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.3f + dec2 * 0.25f + _param.hf_amt * 0.4499f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<6> (comb_fb, lfo1, gains[0]);
    rev.run<7, 8> (comb_fb, flo, glo, fhi, ghi);
    rev.push<6> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<9> (comb_fb.to_const(), overwrite, r, l);

    rev.fetch_block<10> (comb_fb, lfo2, gains[1]);
    rev.run<11, 12> (comb_fb, flo, glo, fhi, ghi);
    rev.push<10> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<13> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<14> (comb_fb, lfo3, gains[2]);
    rev.run<15, 16> (comb_fb, flo, glo, fhi, ghi);
    rev.push<14> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<17> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<18> (comb_fb, lfo4, gains[3]);
    rev.run<19, 20> (comb_fb, flo, glo, fhi, ghi);
    rev.push<18> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<21> (comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    rev.run<22> (l);
    rev.run<23> (l);
    rev.run<24> (comb_fb, in.to_const(), par.character);
    crossfade (l, comb_fb, eramt);

    rev.run<25> (r);
    rev.run<26> (r);
    rev.run<27> (in, par.character);
    crossfade (r, in, eramt);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_dre2000d (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<dre2000d_type> (_modes);

    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, lfo4, tmp1, tmp2;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);
    rev.run<3> (in);
    rev.run<4> (in);
    rev.run<5> (in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 6, 10, 14, 18> (0.25f + dec2 * 5.f, _srate);
    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.80f + _param.lf_amt * 0.199f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.3f + dec2 * 0.3f + _param.hf_amt * 0.29f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<6> (comb_fb, lfo1, gains[0]);
    rev.run<7, 8> (comb_fb, flo, glo, fhi, ghi);
    rev.push<6> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<9> (comb_fb.to_const(), overwrite, r, l);

    rev.fetch_block<10> (comb_fb, blank, gains[1]);
    rev.run<11, 12> (comb_fb, flo, glo, fhi, ghi);
    rev.push<10> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<13> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<14> (comb_fb, blank, gains[2]);
    rev.run<15, 16> (comb_fb, flo, glo, fhi, ghi);
    rev.push<14> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<17> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<18> (comb_fb, lfo4, gains[3]);
    rev.run<19, 20> (comb_fb, flo, glo, fhi, ghi);
    rev.push<18> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<21> (comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    rev.run<22> (l, lfo2);
    rev.run<23> (l);
    rev.run<24> (comb_fb, in.to_const(), par.character);
    crossfade (l, comb_fb, eramt);

    rev.run<25> (r, lfo3);
    rev.run<26> (r);
    rev.run<27> (in, par.character);
    crossfade (r, in, eramt);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_rev5_l_hall (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<rev5_l_hall_type> (_modes);

    block_arr<T> in_mem, main_mem, sub_mem, l_mem, r_mem, lfo1, lfo2, lfo3,
      lfo4, tmp1, tmp2;
    xspan main {main_mem.data(), io.size()};
    xspan sub {sub_mem.data(), io.size()};
    xspan in {in_mem.data(), io.size()};
    xspan l {l_mem.data(), io.size()};
    xspan r {r_mem.data(), io.size()};

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);
    xspan_memdump (main.data(), in);
    xspan_memdump (sub.data(), in);
    rev.run<3> (main, [&] (auto v, uint i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.5_r + 0.2_r * c); // cg3
      tmp2[i] = (T) (0.35_r + 0.7_r * clamp (c, 0._r, 0.4999_r)); // cg2
      return v * 0.25_r;
    });
    rev.run<4> (main, blank, tmp1);
    rev.run<5> (main, blank, tmp2);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 6, 10, 14, 18, 22> (1.f + dec2 * 9.f, _srate);
    T flo = load_float<T> (0.9f + _param.lf_amt * _param.lf_amt * 0.05f);
    T glo = load_float<T> (0.9f + _param.lf_amt * 0.09f);
    T fhi = load_float<T> (0.82f - _param.hf_amt * _param.hf_amt * 0.4f);
    T ghi = load_float<T> (0.6f + _param.hf_amt * 0.35f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<6> (comb_fb, lfo1, gains[0]);
    rev.run<7, 8> (comb_fb, flo, glo, fhi, (T) (ghi * 0.85_r));
    rev.push<6> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<9> (comb_fb.to_const(), overwrite, l, r);

    rev.fetch_block<10> (comb_fb, lfo3, gains[1]);
    rev.run<11, 12> (comb_fb, flo, glo, fhi, (T) (ghi * 0.9_r));
    rev.push<10> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<13> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<14> (comb_fb, blank, gains[2]);
    rev.run<15, 16> (comb_fb, flo, glo, fhi, (T) (ghi * 0.95_r));
    rev.push<14> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<17> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<18> (comb_fb, blank, gains[3]);
    rev.run<19, 20> (comb_fb, flo, glo, fhi, (T) (ghi * 0.97_r));
    rev.push<18> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<21> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<22> (comb_fb, blank, gains[4]);
    rev.run<23, 24> (comb_fb, flo, glo, fhi, ghi);
    rev.push<22> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<25> (comb_fb.to_const(), add_to, l, r);

    rev.run<26> (sub, par.character);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.5_r + 0.2_r * c); // cg3
      // tmp2 still not overwritten...
      // tmp2[i] = (T) (0.35_r + 0.999_r * clamp (c, 0._r, 0.333_r)); // cg2
    }
    rev.run<27> (main, tmp2);
    rev.run<28> (main, tmp1);
    gains
      = rev.get_gain_for_rt60<T, 29, 33, 37, 41, 45> (1.f + dec2 * 5.f, _srate);
    // sub  reverb
    rev.fetch_block<29> (comb_fb, lfo2, gains[0]);
    rev.run<30, 31> (comb_fb, flo, glo, fhi, (T) (ghi * 0.85_r));
    rev.push<29> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<32> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<33> (comb_fb, lfo4, gains[1]);
    rev.run<34, 35> (comb_fb, flo, glo, fhi, (T) (ghi * 0.9_r));
    rev.push<33> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<36> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<37> (comb_fb, blank, gains[2]);
    rev.run<38, 39> (comb_fb, flo, glo, fhi, (T) (ghi * 0.95_r));
    rev.push<37> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<40> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<41> (comb_fb, blank, gains[3]);
    rev.run<42, 43> (comb_fb, flo, glo, fhi, (T) (ghi * 0.97_r));
    rev.push<41> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<44> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<45> (comb_fb, blank, gains[4]);
    rev.run<46, 47> (comb_fb, flo, glo, fhi, ghi);
    rev.push<45> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<48> (comb_fb.to_const(), add_to, l, r);

    rev.run<49> (l, blank, [&] (uint i) {
      return (T) (0.65_r + 0.1_r * lfo1[i]);
    });
    rev.run<50> (l, blank, [&] (uint i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.325_r + 0.999_r * clamp (c, 0._r, 0.333_r)); // cg1
      return (T) (tmp1[i] + 0.05_r * lfo2[i]);
    });
    rev.run<51> (r, blank, [&] (uint i) {
      return (T) (0.65_r + 0.1_r * lfo3[i]);
    });
    rev.run<52> (r, blank, [&] (uint i) {
      return (T) (tmp1[i] + 0.05_r * lfo4[i]);
    });

    rev.run<53> (in.to_const(), add_to, l, r); // ER
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class U, class V>
  T clamp (T v, U min, V max)
  {
    if constexpr (std::is_same_v<T, float>) {
      return std::clamp (v, T {} + min, T {} + max);
    }
    else {
      return fixpt_clamp (v, T {} + min, T {} + max);
    }
  }
  //----------------------------------------------------------------------------
  // 1 selects s2, 0 dst
  template <class T, class U, class V>
  void crossfade (xspan<T> dst, U s2, V ctrl)
  {
    span_visit (dst, [&] (auto spl, uint i) {
      auto ctrl2 = (T) (one<T>() - ctrl[i]);
      dst[i]     = (T) (s2[i] * ctrl[i] + ctrl2 * spl);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T load_float (float v)
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      return fixpt_t::from_float (v);
    }
    else {
      return v;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  static float as_float (T v)
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      return v.to_floatp();
    }
    else {
      return v;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  auto tick_lfo()
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      auto ret = _lfo.tick_sine_fixpt().spec_cast<fixpt_t>().value();
      // At this point "ret" has fixed traits, but different fixed point
      // conversions configured
      return vec_cast<fixpt_t::value_type> (ret);
    }
    else {
      return _lfo.tick_sine();
    }
  }
  //----------------------------------------------------------------------------
  struct mode {
    enum {
      abyss_flt,
      abyss,
      small_space_flt,
      small_space,
      midifex49_flt,
      midifex49,
      midifex50_flt,
      midifex50,
      dre2000a_flt,
      dre2000a,
      dre2000b_flt,
      dre2000b,
      dre2000c_flt,
      dre2000c,
      dre2000d_flt,
      dre2000d,
      rev5_l_hall_flt,
      rev5_l_hall,
#ifdef LOFIVERB_ADD_DEBUG_ALGO
      debug_algo_flt,
      debug_algo,
#endif
    };
  };
  //----------------------------------------------------------------------------
  void update_operation_mode()
  {
    /*
      srate t44k t48k
      7200 49   20
      8400 21   40
      9000 49   16
      10500 21   32
      10800 49   40
      13500 49   32
      14400 49   10
      16800 21   20
      18000 49   08
      20160 35   50
      21000 21   16
      21600 49   20
      22500 49   32
      25200 07   40
      27000 49   16
      28800 49   05
      31500 07   32
      32400 49   40
      33600 21   10
      36000 49   04
      39200 09   60
      39600 49   40
      40320 35   25
      40500 49   32
      42000 21   08
      43200 49   10
      45000 50   16
      46800 52   40
      47040 16   50
      49000 10   49
      49500 55   33
      50400 08   21
      52500 25   35
      54000 60   09
      57600 64   06
      58500 65   39
      58800 04   49
      60480 48   63
      61200 68   51
      63000 10   21
      64800 72   27
      67200 32   07
      67500 75   45
      68400 76   57
    */
    uint srate {};

    switch (_param.mode) {
    case mode::abyss_flt:
    case mode::abyss:
    case mode::small_space_flt:
    case mode::small_space: {
      constexpr auto srates
        = make_array (10500, 16800, 21000, 25200, 31500, 40320, 50400);
      srate = srates[_param.srateid];
    } break;
    case mode::midifex49_flt:
    case mode::midifex49:
    case mode::midifex50_flt:
    case mode::midifex50:
#ifdef LOFIVERB_ADD_DEBUG_ALGO
    case mode::debug_algo_flt:
    case mode::debug_algo:
#endif
    {
      constexpr auto srates
        = make_array (10500, 16800, 21000, 23400, 33600, 42000, 50400);
      srate = srates[_param.srateid];
    } break;
    case mode::dre2000a_flt:
    case mode::dre2000a:
    case mode::dre2000b_flt:
    case mode::dre2000b:
    case mode::dre2000c_flt:
    case mode::dre2000c:
    case mode::dre2000d_flt:
    case mode::dre2000d: {
      constexpr auto srates
        = make_array (16800, 21000, 25200, 32400, 40320, 57600, 57600);
      srate = srates[_param.srateid];
    } break;
    case mode::rev5_l_hall_flt:
    case mode::rev5_l_hall: {
      constexpr auto srates
        = make_array (25200, 33600, 40320, 44100, 45000, 52500, 57600);
      srate = srates[_param.srateid];
    } break;

    default:
      return;
    }
    if (_srate != srate) {
      _srate = srate;
      update_internal_srate (srate, (srate / 20) * 9);
    }

    switch (_param.mode) {
    case mode::abyss_flt:
    case mode::abyss: {
      auto& rev = _modes.emplace<abyss_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    } break;
    case mode::small_space_flt:
    case mode::small_space: {
      auto& rev = _modes.emplace<small_space_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    } break;
    case mode::midifex49_flt:
    case mode::midifex49: {
      auto& rev = _modes.emplace<midifex49_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    } break;
    case mode::midifex50_flt:
    case mode::midifex50: {
      auto& rev = _modes.emplace<midifex50_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    } break;
    case mode::dre2000a_flt:
    case mode::dre2000a: {
      auto& rev = _modes.emplace<dre2000a_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.3333f, 0.6666f, 0.75f});
    } break;
    case mode::dre2000b_flt:
    case mode::dre2000b: {
      auto& rev = _modes.emplace<dre2000b_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.3333f, 0.6666f, 0.75f});
    } break;
    case mode::dre2000c_flt:
    case mode::dre2000c: {
      auto& rev = _modes.emplace<dre2000c_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    } break;
    case mode::dre2000d_flt:
    case mode::dre2000d: {
      auto& rev = _modes.emplace<dre2000d_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    } break;
    case mode::rev5_l_hall_flt:
    case mode::rev5_l_hall: {
      auto& rev = _modes.emplace<rev5_l_hall_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    } break;
#ifdef LOFIVERB_ADD_DEBUG_ALGO
    case mode::debug_algo_flt:
    case mode::debug_algo: {
      auto& rev = _modes.emplace<debug_algo_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    } break;
#endif
    default:
      return;
    }

    _n_processed_samples = 0; // trigger the control block on first sample
    xspan_memset (_mem_reverb, 0);
    update_mod();
  }
  //----------------------------------------------------------------------------
  void update_internal_srate (uint srate, uint fc)
  {
    auto daw_srate = _pc->get_sample_rate();
    auto target_fc = std::min (fc, (daw_srate / 20) * 9);
    _resampler.reset (
      srate,
      daw_srate,
      fc,
      target_fc,
      32,
      16,
      210,
      true,
      max_block_size,
      6 * 1024);

    auto  state   = _pc->get_play_state();
    float beat_hz = 120.f;
    if (state.is_valid) {
      // playhead may not exist, eg. when scanning the plugin
      beat_hz = state.bpm * (1.f / 60.f);
    }
    _1_4beat_spls = (0.25f / beat_hz) * srate;
    _t_spl        = 1.f / srate;

    _ducker.reset();
    _lfo.reset();
    _param_smooth.reset_srate (_t_spl);

    // resize memory
    uint predelay_spls = std::ceil (_1_4beat_spls * max_predelay_qb) * 2;
    uint rev_spls      = 0;
    mp11::mp_for_each<decltype (_modes)> ([&] (auto mode) {
      rev_spls = std::max (mode.get_required_size(), rev_spls);
    });
    uint predelay_flt = predelay_spls * (sizeof (float) / sizeof _mem[0]);

    _mem.clear();
    _mem.resize (predelay_spls + predelay_flt + rev_spls);

    _mem_reverb = xspan {_mem};
    // reminder, hack, _mem_reverb is span<s16> and here it is casting to
    // float. Do keep the predelay (float) first in the memory block to avoid
    // alignment issues.
    _predelay.reset (_mem_reverb.cut_head (predelay_flt).cast<float>(), 2);

    update_mod();
    update_ducker();
  }
  //----------------------------------------------------------------------------
  void update_ducker()
  {
    auto v = _param.ducking_speed * 0.01f;
    _ducker.set_speed (vec_set<2> (v * v), _t_spl);
  }
  //----------------------------------------------------------------------------
  void update_mod()
  {
    // reminder, the phase relations are set on "void set (mode_tag, int v)"
    auto mod = _param_smooth.target().mod;
    switch (_param.mode) {
    case mode::abyss:
    case mode::abyss_flt: {
      auto f_er   = 0.3f + mod * 0.3f;
      auto f_late = 0.1f + mod * 1.2f;
      _lfo.set_freq (f32_x4 {f_er, f_er, f_late, f_late}, _t_spl);
    } break;
    case mode::small_space:
    case mode::small_space_flt: {
      auto f1 = 0.1f + mod * 0.2f;
      _lfo.set_freq (f32_x4 {f1, f1, f1, f1}, _t_spl);
    } break;
    case mode::midifex49_flt:
    case mode::midifex49: {
      auto f = 0.2f + mod * 0.2f;
      _lfo.set_freq (f32_x4 {f, f, f, f}, _t_spl);
    } break;
    case mode::midifex50_flt:
    case mode::midifex50: {
      auto f = 0.3f + mod * 0.1f;
      _lfo.set_freq (f32_x4 {f, f, f, f}, _t_spl);
    } break;
    case mode::dre2000a_flt:
    case mode::dre2000a: {
      auto f = 0.25f + mod * 0.75f;
      _lfo.set_freq (f32_x4 {f, f, f, f}, _t_spl);
    } break;
    case mode::dre2000b_flt:
    case mode::dre2000b: {
      auto f = 0.25f + mod * 0.15f;
      _lfo.set_freq (f32_x4 {f, f, f, f}, _t_spl);
    } break;
    case mode::dre2000c_flt:
    case mode::dre2000c: {
      auto f = 0.25f + mod * 0.75f;
      _lfo.set_freq (f32_x4 {f, f, f, f}, _t_spl);
    } break;
    case mode::dre2000d_flt:
    case mode::dre2000d: {
      auto f = 0.25f + mod * 0.35f;
      _lfo.set_freq (f32_x4 {f, f, f, f}, _t_spl);
    } break;
    case mode::rev5_l_hall_flt:
    case mode::rev5_l_hall: {
      auto f = 0.25f + mod * 0.35f;
      _lfo.set_freq (f32_x4 {f, f, f, f}, _t_spl);
    } break;
    default:
      break;
    }
  }
  //----------------------------------------------------------------------------
  // just a convenience function for iterating block loops while not bloating
  // the code with more loop unroll hints than necessary
  template <class T, class F>
  void span_visit (xspan<T> block, F&& visitor)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block.size(); ++i) {
      visitor (block[i], i);
    }
  }
  //----------------------------------------------------------------------------
  // just a convenience function for iterating block loops while not bloating
  // the code with more loop unroll hints than necessary
  template <class T, class U>
  void span_mul (xspan<T> block, U val)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block.size(); ++i) {
      block[i] = (T) (block[i] * val);
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_max_mode_n_elems()
  {
    uint max_spls {0};
    mp11::mp_for_each<decltype (_modes)> ([&] (auto mode) {
      max_spls = std::max (decltype (mode)::get_required_size(), max_spls);
    });
    return max_spls;
  }
  //----------------------------------------------------------------------------
  template <class T>
  using block_arr = std::array<T, max_block_size>;
  template <class T>
  using fb_block_arr = std::array<T, max_block_size + 1>;
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    u32   mode;
    u32   srateid;
    float gain; // a parameter
    float predelay;
    float lf_amt;
    float hf_amt;
    float ducking_threshold;
    float ducking_speed;
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float mod;
    float stereo;
    float decay;
    float character;
    // dry, wet, ducker/gate
  };
  //----------------------------------------------------------------------------
  struct loop_parameters {
    std::array<fixpt_t, max_block_size> decay;
    std::array<fixpt_t, max_block_size> character;
    std::array<fixpt_t, max_block_size> mod;
    std::array<float, max_block_size>   stereo;
  };
  struct loop_parameters_flt {
    std::array<float, max_block_size> decay;
    std::array<float, max_block_size> character;
    std::array<float, max_block_size> mod;
    std::array<float, max_block_size> stereo;
  };
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  block_resampler<float, 2>             _resampler {};
  static_delay_line<float, true, false> _predelay;

  lfo<4> _lfo;

  using abyss_type
    = detail::lofiverb::engine<detail::lofiverb::abyss_spec, max_block_size>;
  using small_space_type = detail::lofiverb::
    engine<detail::lofiverb::small_space_spec, max_block_size>;
  using midifex49_type = detail::lofiverb::
    engine<detail::lofiverb::midifex49_spec, max_block_size>;
  using midifex50_type = detail::lofiverb::
    engine<detail::lofiverb::midifex50_spec, max_block_size>;
  using dre2000a_type
    = detail::lofiverb::engine<detail::lofiverb::dre2000a_spec, max_block_size>;
  using dre2000b_type
    = detail::lofiverb::engine<detail::lofiverb::dre2000b_spec, max_block_size>;
  using dre2000c_type
    = detail::lofiverb::engine<detail::lofiverb::dre2000c_spec, max_block_size>;
  using dre2000d_type
    = detail::lofiverb::engine<detail::lofiverb::dre2000d_spec, max_block_size>;
  using rev5_l_hall_type = detail::lofiverb::
    engine<detail::lofiverb::rev5_l_hall_spec, max_block_size>;
#ifndef LOFIVERB_ADD_DEBUG_ALGO
  using modes_type = std::variant<
    abyss_type,
    small_space_type,
    midifex49_type,
    midifex50_type,
    dre2000a_type,
    dre2000b_type,
    dre2000c_type,
    dre2000d_type,
    rev5_l_hall_type>;
#else
  using debug_algo_type = detail::lofiverb::
    engine<detail::lofiverb::debug_algo_spec, max_block_size>;

  using modes_type = std::variant<
    debug_algo_type,
    abyss_type,
    small_space_type,
    midifex49_type,
    midifex50_type,
    dre2000a_type,
    dre2000b_type,
    dre2000c_type,
    dre2000d_type,
    rev5_l_hall_type>;
#endif

  modes_type     _modes;
  ducker<f32_x2> _ducker;

  uint  _n_processed_samples;
  float _1_4beat_spls;
  float _t_spl;
  uint  _srate;

  plugin_context* _pc;

  std::vector<s16, overaligned_allocator<s16, 16>> _mem;
  xspan<s16>                                       _mem_reverb;
};
//------------------------------------------------------------------------------
} // namespace artv
