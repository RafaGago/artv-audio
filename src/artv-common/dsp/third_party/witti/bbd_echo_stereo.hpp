#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls

#include "artv-common/dsp/own/classes/add_ducker.hpp"
#include "artv-common/dsp/third_party/witti/bbd_echo.hpp"
#include <algorithm>

namespace artv { namespace witti {

// I liked this one, so I made a low-effort/lazy stereo version. Unfortunately
// The parameter count was already high and I only had 2 sliders left, so one
// is for the synced L/R and the other one became a decalibration parameter to
// create stereo differences.
//------------------------------------------------------------------------------
struct bbd_echo_stereo : private add_ducker<f64_x2> {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::delay;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    if (v == _params.feedback) {
      return;
    }
    _params.feedback = v;
    update();
  }
  static constexpr auto get_parameter (feedback_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::feedback_tag {});
  }
  //----------------------------------------------------------------------------
  struct lfo_depth_tag {};
  void set (lfo_depth_tag, float v)
  {
    if (v == _params.lfo_depth) {
      return;
    }
    _params.lfo_depth = v;
    update();
  }

  static constexpr auto get_parameter (lfo_depth_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::lfo_depth_tag {});
  }
  //----------------------------------------------------------------------------
  struct lfo_speed_tag {};
  void set (lfo_speed_tag, float v)
  {
    if (v == _params.lfo_speed) {
      return;
    }
    _params.lfo_speed = v;
    update();
  }

  static constexpr auto get_parameter (lfo_speed_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::lfo_speed_tag {});
  }
  //----------------------------------------------------------------------------
  struct stages_tag {};
  void set (stages_tag, int v)
  {
    if (v == _params.stages) {
      return;
    }
    _params.stages = v;
    update();
  }
  static constexpr auto get_parameter (stages_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::stages_tag {});
  }
  //----------------------------------------------------------------------------
  struct delay_tag {};
  void set (delay_tag, float v)
  {
    if (v == _params.delay) {
      return;
    }
    _params.delay = v;
    update();
  }
  static constexpr auto get_parameter (delay_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::delay_tag {});
  }
  //----------------------------------------------------------------------------
  struct delay_sync_l_tag {};
  void set (delay_sync_l_tag, float v)
  {
    if (v == _params.delay_sync_l) {
      return;
    }
    _params.delay_sync_l = v;
    update();
  }

  static constexpr auto get_parameter (delay_sync_l_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::delay_sync_tag {});
  }
  //----------------------------------------------------------------------------
  struct delay_sync_r_tag {};
  void set (delay_sync_r_tag, float v)
  {
    if (v == _params.delay_sync_r) {
      return;
    }
    _params.delay_sync_r = v;
    update();
  }

  static constexpr auto get_parameter (delay_sync_r_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::delay_sync_tag {});
  }
  //----------------------------------------------------------------------------
  struct hp_filter_tag {};
  void set (hp_filter_tag, float v)
  {
    if (v == _params.hp_filter) {
      return;
    }
    _params.hp_filter = v;
    update();
  }
  static constexpr auto get_parameter (hp_filter_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::hp_filter_tag {});
  }
  //----------------------------------------------------------------------------
  struct hp_res_tag {};
  void set (hp_res_tag, float v)
  {
    if (v == _params.hp_res) {
      return;
    }
    _params.hp_res = v;
    update();
  }
  static constexpr auto get_parameter (hp_res_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::hp_res_tag {});
  }
  //----------------------------------------------------------------------------
  struct lp_filter_tag {};
  void set (lp_filter_tag, float v)
  {
    if (v == _params.lp_filter) {
      return;
    }
    _params.lp_filter = v;
    update();
  }
  static constexpr auto get_parameter (lp_filter_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::lp_filter_tag {});
  }
  //----------------------------------------------------------------------------
  struct lp_res_tag {};
  void set (lp_res_tag, float v)
  {
    if (v == _params.lp_res) {
      return;
    }
    _params.lp_res = v;
    update();
  }
  static constexpr auto get_parameter (lp_res_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::lp_res_tag {});
  }
  //----------------------------------------------------------------------------
  struct clock_offset_tag {};
  void set (clock_offset_tag, float v)
  {
    if (v == _params.clock_offset) {
      return;
    }
    _params.clock_offset = v;
    update();
  }
  static constexpr auto get_parameter (clock_offset_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::clock_offset_tag {});
  }
  //----------------------------------------------------------------------------
  struct clock_scale_tag {};
  void set (clock_scale_tag, float v)
  {
    if (v == _params.clock_scale) {
      return;
    }
    _params.clock_scale = v;
    update();
  }
  static constexpr auto get_parameter (clock_scale_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::clock_scale_tag {});
  }
  //----------------------------------------------------------------------------
  struct clock_curve_tag {};
  void set (clock_curve_tag, float v)
  {
    if (v == _params.clock_curve) {
      return;
    }
    _params.clock_curve = v;
    update();
  }
  static constexpr auto get_parameter (clock_curve_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::clock_curve_tag {});
  }
  //----------------------------------------------------------------------------
  struct hiss_tag {};
  void set (hiss_tag, float v)
  {
    if (v == _params.hiss) {
      return;
    }
    _params.hiss = v;
    update();
  }
  static constexpr auto get_parameter (hiss_tag)
  {
    return bbd_echo::get_parameter (bbd_echo::hiss_tag {});
  }
  //----------------------------------------------------------------------------
  struct decalibration_tag {};
  void set (decalibration_tag, float v)
  {
    v = v * 0.01;
    if (v == _params.decalibration) {
      return;
    }
    _params.decalibration = v;
    update();
  }
  static constexpr auto get_parameter (decalibration_tag)
  {
    return float_param ("", -100.0, 100.0, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  using add_ducker::get_parameter;
  using add_ducker::set;
  using ducking_speed_tag     = add_ducker::ducking_speed_tag;
  using ducking_threshold_tag = add_ducker::ducking_threshold_tag;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    feedback_tag,
    lfo_depth_tag,
    lfo_speed_tag,
    stages_tag,
    delay_tag,
    delay_sync_l_tag,
    delay_sync_r_tag,
    hp_filter_tag,
    hp_res_tag,
    lp_filter_tag,
    lp_res_tag,
    clock_offset_tag,
    clock_scale_tag,
    clock_curve_tag,
    hiss_tag,
    ducking_speed_tag,
    ducking_threshold_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    add_ducker::reset (pc.get_sample_rate());
    _params.feedback
      = bbd_echo::get_parameter (bbd_echo::feedback_tag {}).defaultv;
    _params.lfo_depth
      = bbd_echo::get_parameter (bbd_echo::lfo_depth_tag {}).defaultv;
    _params.lfo_speed
      = bbd_echo::get_parameter (bbd_echo::lfo_speed_tag {}).defaultv;
    _params.delay = bbd_echo::get_parameter (bbd_echo::delay_tag {}).defaultv;
    _params.hp_filter
      = bbd_echo::get_parameter (bbd_echo::hp_filter_tag {}).defaultv;
    _params.hp_filter
      = bbd_echo::get_parameter (bbd_echo::lp_filter_tag {}).defaultv;
    _params.hp_res = bbd_echo::get_parameter (bbd_echo::hp_res_tag {}).defaultv;
    _params.hp_res = bbd_echo::get_parameter (bbd_echo::lp_res_tag {}).defaultv;
    _params.clock_offset
      = bbd_echo::get_parameter (bbd_echo::clock_offset_tag {}).defaultv;
    _params.clock_scale
      = bbd_echo::get_parameter (bbd_echo::clock_scale_tag {}).defaultv;
    _params.clock_curve
      = bbd_echo::get_parameter (bbd_echo::clock_curve_tag {}).defaultv;
    _params.hiss = bbd_echo::get_parameter (bbd_echo::hiss_tag {}).defaultv;

    _params.delay_sync_l = _params.delay_sync_r
      = bbd_echo::get_parameter (bbd_echo::delay_sync_tag {}).defaultv;
    _params.stages = bbd_echo::get_parameter (bbd_echo::stages_tag {}).defaultv;

    _l.reset (pc);
    _r.reset (pc);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    add_ducker::process (
      outs,
      ins,
      samples,
      [=] (xspan<T*> outs_fw, xspan<T const*> ins_fw, uint samples_fw) {
        this->process_intern (outs_fw, ins_fw, samples_fw);
      });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  void process_intern (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint i = 0; i < (uint) bus_type; ++i) {
      if (unlikely (ins[i] != outs[i])) {
        memcpy (outs[i], ins[i], block_samples * sizeof outs[i][0]);
      }
    }
    _l.process (make_array (outs[0]), block_samples);
    _r.process (make_array (outs[1]), block_samples);
  }
  //----------------------------------------------------------------------------
  void update()
  {
    float decal = _params.decalibration * 0.1f;

    float decal_ha = 1.f + decal;
    float decal_hb = 1.f - decal;
    float decal_ma = 1.f + (decal * 0.2);
    float decal_mb = 1.f + -(decal * 0.1);
    float decal_la = 1.f + (decal * 0.01);
    float decal_lb = 1.f - (decal * 0.02);

    set<bbd_echo::feedback_tag> (_params.feedback, decal_ha, decal_hb);
    set<bbd_echo::lfo_depth_tag> (_params.lfo_depth, decal_hb, decal_ha);
    set<bbd_echo::lfo_speed_tag> (_params.lfo_speed, decal_la, decal_lb);
    set<bbd_echo::delay_tag> (_params.delay, decal_lb, decal_la);
    set<bbd_echo::hp_filter_tag> (_params.hp_filter, decal_mb, decal_ma);
    set<bbd_echo::lp_filter_tag> (_params.lp_filter, decal_mb, decal_ma);
    set<bbd_echo::hp_res_tag> (_params.hp_res, decal_ha, decal_hb);
    set<bbd_echo::lp_res_tag> (_params.lp_res, decal_hb, decal_ha);
    set<bbd_echo::clock_offset_tag> (_params.clock_offset, decal_hb, decal_ha);
    set<bbd_echo::clock_scale_tag> (_params.clock_scale, decal_ha, decal_hb);
    set<bbd_echo::clock_curve_tag> (_params.clock_curve, decal_hb, decal_ha);
    set<bbd_echo::hiss_tag> (_params.hiss, decal_hb, decal_ha);

    _l.set (bbd_echo::delay_sync_tag {}, _params.delay_sync_l);
    _r.set (bbd_echo::delay_sync_tag {}, _params.delay_sync_r);
    _l.set (bbd_echo::stages_tag {}, _params.stages);
    _r.set (bbd_echo::stages_tag {}, _params.stages);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void set (float v, float decal_l, float decal_r)
  {
    constexpr auto min = bbd_echo::get_parameter (T {}).min;
    constexpr auto max = bbd_echo::get_parameter (T {}).max;
    _l.set (T {}, std::clamp (v * decal_l, min, max));
    _r.set (T {}, std::clamp (v * decal_r, min, max));
  }
  //----------------------------------------------------------------------------
  struct all_params {
    float feedback;
    float lfo_depth;
    float lfo_speed;
    int   stages;
    float delay;
    float delay_sync_l;
    float delay_sync_r;
    float hp_filter;
    float hp_res;
    float lp_filter;
    float lp_res;
    float clock_offset;
    float clock_scale;
    float clock_curve;
    float hiss;
    float decalibration;
  };
  //----------------------------------------------------------------------------
  all_params _params;
  bbd_echo   _l;
  bbd_echo   _r;
};
//------------------------------------------------------------------------------
}} // namespace artv::witti
