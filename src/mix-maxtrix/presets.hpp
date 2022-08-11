#pragma once

namespace artv {

struct preset {
  char const* name;
  char const* xml;
};

// clang-format off
static constexpr preset presets[] = {
    {.name = "Routing: Bus 1 Volume/Pan/Width control",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
</params>
  )END"},
  {.name = "Routing: Passhtrough",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="in_selection_02" value="2.0" />
    <PARAM id="in_selection_03" value="4.0" />
    <PARAM id="in_selection_04" value="8.0" />
    <PARAM id="in_selection_05" value="16.0" />
    <PARAM id="in_selection_06" value="32.0" />
    <PARAM id="in_selection_07" value="64.0" />
    <PARAM id="in_selection_08" value="128.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="out_selection_02" value="2.0" />
    <PARAM id="out_selection_03" value="4.0" />
    <PARAM id="out_selection_04" value="8.0" />
    <PARAM id="out_selection_05" value="16.0" />
    <PARAM id="out_selection_06" value="32.0" />
    <PARAM id="out_selection_07" value="64.0" />
    <PARAM id="out_selection_08" value="128.0" />
</params>
  )END"},
  {.name = "Routing: Bus 1 Joiner",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="255.0" />
    <PARAM id="out_selection_01" value="1.0" />
</params>
  )END"},
  {.name = "Routing: Bus 1 Forwarder",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="255.0" />
</params>
  )END"},
  {.name = "Routing: Swap Bus1 with Bus2",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="in_selection_02" value="2.0" />
    <PARAM id="in_selection_03" value="4.0" />
    <PARAM id="in_selection_04" value="8.0" />
    <PARAM id="in_selection_05" value="16.0" />
    <PARAM id="in_selection_06" value="32.0" />
    <PARAM id="in_selection_07" value="64.0" />
    <PARAM id="in_selection_08" value="128.0" />
    <PARAM id="out_selection_01" value="2.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="out_selection_03" value="4.0" />
    <PARAM id="out_selection_04" value="8.0" />
    <PARAM id="out_selection_05" value="16.0" />
    <PARAM id="out_selection_06" value="32.0" />
    <PARAM id="out_selection_07" value="64.0" />
    <PARAM id="out_selection_08" value="128.0" />
</params>
     )END"},
  {.name = "Mix: 4-band Linkwitz-Riley crossover",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="lr_crossv_band1_diff" />
    <PARAM id="lr_crossv_band1_frequency" value="51.08850860595703" />
    <PARAM id="lr_crossv_band1_mode" />
    <PARAM id="lr_crossv_band1_out" value="0.0" />
    <PARAM id="lr_crossv_band2_diff" />
    <PARAM id="lr_crossv_band2_frequency" value="73.13029479980469" />
    <PARAM id="lr_crossv_band2_mode" />
    <PARAM id="lr_crossv_band2_out" value="1.0" />
    <PARAM id="lr_crossv_band3_diff" />
    <PARAM id="lr_crossv_band3_frequency" value="105.9929809570312" />
    <PARAM id="lr_crossv_band3_mode" />
    <PARAM id="lr_crossv_band3_out" value="2.0" />
    <PARAM id="fx_type_01" value="42.0" />
    <PARAM id="global_volume" />
    <PARAM id="mixer_sends" />
    <PARAM id="out_selection_02" value="2.0" />
    <PARAM id="out_selection_03" value="4.0" />
    <PARAM id="out_selection_04" value="8.0" />
    <PARAM id="out_selection_05" value="0.0" />
    <PARAM id="out_selection_07" value="0.0" />
    <PARAM id="routing" value="0.0" />
    <gui_data  bus1_text="Lows" bus2_text="Low Mids" bus3_text="Mids" bus4_text="Highs" />
</params>
     )END"},
  {.name = "Mix: Console 6 Companding",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="255.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="consoles_type_01" value="2.0" />
    <PARAM id="fx_type_01" value="1.0" />
</params>
     )END"},
  {.name = "Mix: Channel Strip",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="0.0" />
    <PARAM id="fx_type_01" value="22.0" />
    <PARAM id="fx_type_02" value="16.0" />
    <PARAM id="fx_type_03" value="39.0" />
    <PARAM id="fx_type_04" value="11.0" />
    <PARAM id="fx_type_05" value="20.0" />
    <PARAM id="global_volume" />
    <PARAM id="mixer_sends" value="15.0" />
    <PARAM id="out_selection_05" value="1.0" />
</params>
     )END"},
  {.name = "Mix: EQ: Left, Right, Mid, Side, Full.",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="fx_type_01" value="23.0" />
    <PARAM id="fx_type_02" value="23.0" />
    <PARAM id="fx_type_03" value="23.0" />
    <PARAM id="fx_type_04" value="23.0" />
    <PARAM id="fx_type_05" value="23.0" />
    <PARAM id="fx_type_06" value="0.0" />
    <PARAM id="in_selection_02" value="1.0" />
    <PARAM id="in_selection_03" value="1.0" />
    <PARAM id="in_selection_04" value="1.0" />
    <PARAM id="in_selection_05" value="1.0" />
    <PARAM id="in_selection_06" value="0.0" />
    <PARAM id="in_selection_07" value="1.0" />
    <PARAM id="in_selection_08" value="0.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="out_selection_03" value="1.0" />
    <PARAM id="out_selection_04" value="1.0" />
    <PARAM id="out_selection_05" value="1.0" />
    <PARAM id="out_selection_06" value="0.0" />
    <PARAM id="out_selection_07" value="1.0" />
    <PARAM id="out_selection_08" value="0.0" />
    <PARAM id="routing" value="2.0" />
    <PARAM id="wet_balance_03" value="-100.0" />
    <PARAM id="wet_balance_04" value="100.0" />
    <PARAM id="wet_pan_01" value="-100.0" />
    <PARAM id="wet_pan_02" value="100.0" />
    <gui_data bus1_text="Left" bus2_text="Right" bus3_text="Mid" bus4_text="Side" bus5_text="All" bus7_text="Passhtrough" />
</params>
     )END"},
  {.name = "Mix: Transient Splitter",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1002000">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="fx_type_01" value="45.0" />
    <PARAM id="mixer_sends" value="128.0" />
    <PARAM id="out_selection_02" value="2.0" />
    <PARAM id="transient_gate_decay_01" value="47.40000152587891" />
    <PARAM id="transient_gate_decay_shape_01" value="0.0" />
    <PARAM id="transient_gate_detect_hipass_01" value="39.64512634277344" />
    <PARAM id="transient_gate_detect_recovery_01" value="29.39999961853027" />
    <PARAM id="transient_gate_detect_shape_01" value="0.0" />
    <gui_data bus1_text="Transient" bus2_text="Tonal" bus3_text="" bus4_text="" bus5_text="" bus6_text="" bus7_text="" bus8_text="" />
</params>
    )END"},
  {.name = "Mix: Mojo and Tails",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1004000">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="0.0" />
    <PARAM id="eq4x_band1_freq_01" value="41.90177917480469" />
    <PARAM id="eq4x_band1_freq_03" value="110.6799926757812" />
    <PARAM id="eq4x_band1_gain_01" value="2.25" />
    <PARAM id="eq4x_band1_gain_02" value="0.0" />
    <PARAM id="eq4x_band1_gain_03" value="9.0" />
    <PARAM id="eq4x_band1_q_01" value="1.570600032806396" />
    <PARAM id="eq4x_band1_q_02" value="0.4999999701976776" />
    <PARAM id="eq4x_band1_q_03" value="0.100100003182888" />
    <PARAM id="eq4x_band1_type_01" value="10.0" />
    <PARAM id="eq4x_band1_type_02" value="0.0" />
    <PARAM id="eq4x_band1_type_03" value="10.0" />
    <PARAM id="eq4x_band2_freq_02" value="69.0" />
    <PARAM id="fx_type_01" value="23.0" />
    <PARAM id="fx_type_02" value="41.0" />
    <PARAM id="fx_type_03" value="23.0" />
    <PARAM id="fx_type_04" value="45.0" />
    <PARAM id="in_selection_03" value="1.0" />
    <PARAM id="in_selection_08" value="1.0" />
    <PARAM id="mixer_sends" value="5.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="out_selection_04" value="1.0" />
    <PARAM id="out_selection_08" value="1.0" />
    <PARAM id="wet_balance_02" value="-70.40000152587891" />
    <PARAM id="wet_balance_04" value="60.0" />
    <PARAM id="transient_gate_decay_04" value="61.40000152587891" />
    <PARAM id="transient_gate_detect_channels_04" value="1.0" />
    <PARAM id="transient_gate_detect_hipass_04" value="57.50381469726562" />
    <PARAM id="transient_gate_detect_recovery_04" value="11.40000057220459" />
    <PARAM id="transient_gate_detect_shape_04" value="2.0" />
    <PARAM id="transient_gate_output_04" value="1.0" />
    <PARAM id="waveshaper_drive_02" value="17.25" />
    <PARAM id="waveshaper_drive_balance_02" value="-1.75" />
    <PARAM id="waveshaper_envfollow_attack_02" value="9.960000038146973" />
    <PARAM id="waveshaper_envfollow_mode_02" value="2.0" />
    <PARAM id="waveshaper_envfollow_release_02" value="83.51999664306641" />
    <PARAM id="waveshaper_envfollow_sensitivity_02" value="20.0" />
    <PARAM id="waveshaper_envfollow_to_dc_02" value="55.19999694824219" />
    <PARAM id="waveshaper_envfollow_to_drive_02" value="-70.40000152587891" />
    <PARAM id="waveshaper_feedback_02" value="8.300003051757812" />
    <PARAM id="waveshaper_mode_02" value="1.0" />
    <PARAM id="waveshaper_type_02" value="2.0" />
    <gui_data bus1_text="Lows" bus2_text="Low mojo" bus3_text="Highs" bus4_text="Extra tails" bus5_text="" bus6_text="" bus7_text="" bus8_text="Mute me!" />
</params>
     )END"},
  {.name = "Mix: Long Sparse Reverb (metals, sax)",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1005000">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="fx_type_01" value="50.0" />
    <PARAM id="in_selection_02" value="1.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="reverb_damp_factor_01" value="40.0" />
    <PARAM id="reverb_damp_freq_01" value="-56.79999923706055" />
    <PARAM id="reverb_diff_in_01" value="0.0" />
    <PARAM id="reverb_diff_out_01" value="0.0" />
    <PARAM id="reverb_early_2_late_bal_01" value="47.5" />
    <PARAM id="reverb_early_gain_01" value="-15.0" />
    <PARAM id="reverb_early_size_01" value="-78.40000152587891" />
    <PARAM id="reverb_gap_01" value="0.0" />
    <PARAM id="reverb_hp_freq_01" value="48.59999847412109" />
    <PARAM id="reverb_l_sparseness_01" value="-55.26000213623047" />
    <PARAM id="reverb_late_gain_01" value="-9.0" />
    <PARAM id="reverb_lf_time_factor_01" value="24.0" />
    <PARAM id="reverb_lr_sparseness_01" value="-73.16999816894531" />
    <PARAM id="reverb_mod_depth_01" value="63.34000015258789" />
    <PARAM id="reverb_mod_freq_01" value="0.3920000195503235" />
    <PARAM id="reverb_mod_wave_01" value="5.0" />
    <PARAM id="reverb_r_sparseness_01" value="-65.69000244140625" />
    <PARAM id="reverb_size_01" value="72.80000305175781" />
    <PARAM id="reverb_decay_01" value="12133.0" />
    <gui_data bus1_text="Wet" bus2_text="Dry" bus3_text="" bus4_text="" bus5_text="" bus6_text="" bus7_text="" bus8_text="" />
</params>
   )END"},
  {.name = "Mix: Long Reverb",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1004000">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="fx_type_01" value="39.0" />
    <PARAM id="fx_type_02" value="50.0" />
    <PARAM id="in_selection_03" value="1.0" />
    <PARAM id="mixer_sends" value="1.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="out_selection_03" value="1.0" />
    <PARAM id="reverb_damp_factor_02" value="30.0" />
    <PARAM id="reverb_damp_freq_02" value="-4.800003051757812" />
    <PARAM id="reverb_diff_in_02" value="32.0" />
    <PARAM id="reverb_diff_out_02" value="100.0" />
    <PARAM id="reverb_early_gain_02" value="-11.0" />
    <PARAM id="reverb_early_size_02" value="11.74100494384766" />
    <PARAM id="reverb_gap_02" value="2.944000244140625" />
    <PARAM id="reverb_hp_freq_02" value="60.46999740600586" />
    <PARAM id="reverb_l_sparseness_02" value="-2.910003662109375" />
    <PARAM id="reverb_late_gain_02" value="-7.0" />
    <PARAM id="reverb_lf_time_factor_02" value="48.0" />
    <PARAM id="reverb_lr_sparseness_02" value="1.349998474121094" />
    <PARAM id="reverb_mod_depth_02" value="42.39999771118164" />
    <PARAM id="reverb_mod_freq_02" value="22.07600021362305" />
    <PARAM id="reverb_mod_wave_02" value="1.0" />
    <PARAM id="reverb_predelay_02" value="0.4880000352859497" />
    <PARAM id="reverb_r_sparseness_02" value="-0.7900009155273438" />
    <PARAM id="reverb_size_02" value="7.62939453125e-6" />
    <PARAM id="reverb_stereo_02" value="84.79999542236328" />
    <PARAM id="reverb_time_02" value="20000.0" />
    <PARAM id="track_comp_compattack_01" value="1.85099995136261" />
    <PARAM id="track_comp_compfeedbk_01" value="0.0" />
    <PARAM id="track_comp_compratio_01" value="6.200000286102295" />
    <PARAM id="track_comp_comprelease_01" value="137.3099975585938" />
    <PARAM id="track_comp_compthresh_01" value="-34.31999969482422" />
    <PARAM id="track_comp_saturation_01" value="0.0" />
    <PARAM id="fx_mix_02" value="100.0" />
    <PARAM id="reverb_damp_freq_01" value="0.0" />
    <PARAM id="reverb_decay_01" value="1500.0" />
    <PARAM id="reverb_decay_02" value="19017.0" />
    <PARAM id="reverb_early_2_late_bal_01" value="0.0" />
    <PARAM id="reverb_early_2_late_bal_02" value="0.0" />
    <PARAM id="reverb_early_gain_01" value="-6.0" />
    <PARAM id="reverb_early_size_01" value="7.62939453125e-6" />
    <PARAM id="reverb_gap_01" value="3.000000238418579" />
    <PARAM id="reverb_hp_freq_01" value="35.0" />
    <PARAM id="reverb_l_sparseness_01" value="2.339996337890625" />
    <PARAM id="reverb_late_gain_01" value="-6.0" />
    <PARAM id="reverb_lf_time_factor_01" value="0.0" />
    <PARAM id="reverb_lr_sparseness_01" value="1.229995727539062" />
    <PARAM id="reverb_mod_depth_01" value="10.0" />
    <PARAM id="reverb_mod_freq_01" value="50.00000381469727" />
    <PARAM id="reverb_mod_mode_02" value="0.0" />
    <PARAM id="reverb_r_sparseness_01" value="-1.0" />
    <PARAM id="reverb_size_01" value="7.62939453125e-6" />
    <gui_data bus1_text="Pre-Comp" bus2_text="Wet" bus3_text="Dry" bus4_text="-" bus5_text="" bus6_text="" bus7_text="" bus8_text="" />
</params>
     )END"},
  {.name = "FX: Phased Delay",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="in_selection_02" value="1.0" />
    <PARAM id="in_selection_04" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="out_selection_03" value="1.0" />
    <PARAM id="out_selection_04" value="1.0" />
    <PARAM id="bbd_echo_age_02" value="-61.59999847412109" />
    <PARAM id="bbd_echo_clock_curve_02" value="1.039999961853027" />
    <PARAM id="bbd_echo_clock_offset_02" value="15.20100116729736" />
    <PARAM id="bbd_echo_clock_scale_02" value="0.3950000107288361" />
    <PARAM id="bbd_echo_delay_02" value="0.0" />
    <PARAM id="bbd_echo_delay_sync_l_02" value="25.19999885559082" />
    <PARAM id="bbd_echo_delay_sync_r_02" value="46.39999771118164" />
    <PARAM id="bbd_echo_feedback_02" value="31.19999694824219" />
    <PARAM id="bbd_echo_hp_filter_02" value="33.0" />
    <PARAM id="bbd_echo_hp_res_02" value="0.5299999713897705" />
    <PARAM id="bbd_echo_lfo_depth_02" value="0.02999999932944775" />
    <PARAM id="bbd_echo_lfo_speed_02" value="0.07999999821186066" />
    <PARAM id="bbd_echo_lp_filter_02" value="89.0" />
    <PARAM id="bbd_echo_lp_res_02" value="0.6800000071525574" />
    <PARAM id="df_er_high_cut_04" value="115.6610946655273" />
    <PARAM id="df_er_low_cut_04" value="-29.10841751098633" />
    <PARAM id="df_er_size_04" value="13.39999961853027" />
    <PARAM id="df_er_width_04" value="75.19999694824219" />
    <PARAM id="dry_balance_03" value="0.0" />
    <PARAM id="echo_cycles_feedback_ratio_02" value="72.79999542236328" />
    <PARAM id="echo_cycles_feedback_rotation_02" value="166.3000030517578" />
    <PARAM id="echo_cycles_filter_bandwidth_02" value="3.549999952316284" />
    <PARAM id="echo_cycles_filter_db_02" value="6.400000095367432" />
    <PARAM id="echo_cycles_filter_freq_02" value="68.71851348876953" />
    <PARAM id="echo_cycles_input_width_02" value="49.20000076293945" />
    <PARAM id="echo_cycles_output_variation_02" value="23.60000038146973" />
    <PARAM id="echo_cycles_rotation_mode_02" value="1.0" />
    <PARAM id="fx_mix_03" value="73.20000457763672" />
    <PARAM id="fx_type_02" value="4.0" />
    <PARAM id="fx_type_03" value="25.0" />
    <PARAM id="fx_type_04" value="30.0" />
    <PARAM id="mixer_sends" value="2.0" />
    <PARAM id="myphaser_feedback_03" value="-72.0" />
    <PARAM id="myphaser_high_freq_03" value="66.44156646728516" />
    <PARAM id="myphaser_lfo_depth_03" value="53.79999923706055" />
    <PARAM id="myphaser_lfo_rate_03" value="-80.0" />
    <PARAM id="myphaser_lfo_rate_sync_03" value="3.0" />
    <PARAM id="myphaser_lfo_stereo_03" value="0.0" />
    <PARAM id="myphaser_low_freq_03" value="45.76940155029297" />
    <PARAM id="myphaser_q_03" value="0.2470000088214874" />
    <PARAM id="myphaser_stages_03" value="11.0" />
    <PARAM id="myphaser_stages_mode_03" value="4.0" />
    <PARAM id="volume_03" value="-9.900001525878906" />
    <PARAM id="volume_04" value="-12.97000122070312" />
    <PARAM id="wet_balance_02" value="0.0" />
    <PARAM id="wet_balance_03" value="0.0" />
    <PARAM id="wet_pan_03" value="8.799995422363281" />
    <PARAM id="wet_pan_04" value="-1.599998474121094" />
</params>
     )END"},
  {.name = "FX: Filter Rumble",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="out_selection_04" value="1.0" />
    <PARAM id="fx_mix_01" value="100.0" />
    <PARAM id="fx_mix_02" value="100.0" />
    <PARAM id="fx_type_01" value="0.0" />
    <PARAM id="fx_type_02" value="37.0" />
    <PARAM id="fx_type_03" value="28.0" />
    <PARAM id="global_volume" value="-6.0" />
    <PARAM id="mixer_sends" value="7.0" />
    <PARAM id="mute_solo_02" value="0.0" />
    <PARAM id="volume_01" value="-1.080001831054688" />
    <PARAM id="volume_02" value="1.520000457763672" />
    <PARAM id="volume_03" value="3.149997711181641" />
    <PARAM id="wet_pan_01" value="0.0" />
    <PARAM id="wet_pan_02" value="0.0" />
    <PARAM id="wet_pan_03" value="48.0" />
    <PARAM id="channel_modifs_02" value="0.0" />
    <PARAM id="channel_modifs_04" value="0.0" />
    <PARAM id="df_room_decay_04" value="0.9000000357627869" />
    <PARAM id="df_room_diffuse_04" value="19.20000076293945" />
    <PARAM id="df_room_early_04" value="-56.20000076293945" />
    <PARAM id="df_room_early_damp_04" value="119.1286392211914" />
    <PARAM id="df_room_early_send_04" value="-16.29999923706055" />
    <PARAM id="df_room_late_04" value="0.0" />
    <PARAM id="df_room_late_damp_04" value="103.204231262207" />
    <PARAM id="df_room_low_boost_04" value="37.0" />
    <PARAM id="df_room_low_boost_freq_04" value="72.79496002197266" />
    <PARAM id="df_room_low_cut_04" value="55.34995651245117" />
    <PARAM id="df_room_size_04" value="27.0" />
    <PARAM id="df_room_spin_04" value="-64.52828216552734" />
    <PARAM id="df_room_wander_04" value="68.0" />
    <PARAM id="filter2x_band1_drive_02" value="-26.79999923706055" />
    <PARAM id="filter2x_band1_feedback_02" value="-71.19999694824219" />
    <PARAM id="filter2x_band1_freq_02" value="57.51941680908203" />
    <PARAM id="filter2x_band1_gain_02" value="16.20000076293945" />
    <PARAM id="filter2x_band1_reso_02" value="0.6507999897003174" />
    <PARAM id="filter2x_band1_tolerance_02" value="40.80000305175781" />
    <PARAM id="filter2x_band1_topology_02" value="2.0" />
    <PARAM id="filter2x_band1_type_02" value="3.0" />
    <PARAM id="filter2x_band2_drive_02" value="-18.09999847412109" />
    <PARAM id="filter2x_band2_feedback_02" value="-41.59999847412109" />
    <PARAM id="filter2x_band2_freq_02" value="58.47613143920898" />
    <PARAM id="filter2x_band2_gain_02" value="9.0" />
    <PARAM id="filter2x_band2_reso_02" value="0.7960000038146973" />
    <PARAM id="filter2x_band2_tolerance_02" value="-38.39999771118164" />
    <PARAM id="filter2x_band2_topology_02" value="3.0" />
    <PARAM id="filter2x_band2_type_02" value="1.0" />
    <PARAM id="fx_mix_03" value="17.60000038146973" />
    <PARAM id="fx_mix_04" value="8.40000057220459" />
    <PARAM id="fx_type_04" value="33.0" />
    <PARAM id="volume_04" value="4.30999755859375" />
    <PARAM id="wet_balance_02" value="0.0" />
    <PARAM id="wet_pan_04" value="45.59999084472656" />
    <PARAM id="zelittlechorus_depth_03" value="11.90399932861328" />
    <PARAM id="zelittlechorus_f0_03" value="3246.39990234375" />
    <PARAM id="zelittlechorus_f1_03" value="2058.400146484375" />
    <PARAM id="zelittlechorus_fb_03" value="0.3720000088214874" />
    <PARAM id="zelittlechorus_fblp_03" value="0.2899999916553497" />
    <PARAM id="zelittlechorus_fbtype_03" value="1.0" />
    <PARAM id="zelittlechorus_lfoa_03" value="10.94999980926514" />
    <PARAM id="zelittlechorus_lfob_03" value="7.423799991607666" />
    <PARAM id="zelittlechorus_lfof_03" value="6.364999771118164" />
    <PARAM id="zelittlechorus_lfomix_03" value="0.4679999947547913" />
    <PARAM id="zelittlechorus_r_03" value="0.7965999841690063" />
    <PARAM id="zelittlechorus_selfpm_03" value="0.812000036239624" />
</params>
     )END"},
  {.name = "FX: Hard Panned Plate",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="dry_pan_01" value="100.0" />
    <PARAM id="fx_mix_01" value="95.20000457763672" />
    <PARAM id="fx_type_01" value="34.0" />
    <PARAM id="tal_reverb2_decay_01" value="27.60000038146973" />
    <PARAM id="tal_reverb2_highshelf_frequency_01" value="114.6010284423828" />
    <PARAM id="tal_reverb2_highshelf_gain_01" value="9.220000267028809" />
    <PARAM id="tal_reverb2_lowshelf_frequency_01" value="55.46834945678711" />
    <PARAM id="tal_reverb2_lowshelf_gain_01" value="9.210000038146973" />
    <PARAM id="tal_reverb2_peak_frequency_01" value="60.90897750854492" />
    <PARAM id="tal_reverb2_peak_gain_01" value="16.42000007629395" />
    <PARAM id="tal_reverb2_stereo_width_01" value="10.40000057220459" />
    <PARAM id="wet_pan_01" value="-100.0" />
</params>
     )END"},
  {.name = "FX: Vocal Ambience",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_04" value="1.0" />
    <PARAM id="_4x4_high_drive_02" value="72.80000305175781" />
    <PARAM id="_4x4_high_gain_02" value="4.30000114440918" />
    <PARAM id="_4x4_low_mid_freq_02" value="57.91906356811523" />
    <PARAM id="_4x4_mid_drive_02" value="75.20000457763672" />
    <PARAM id="_4x4_mid_high_freq_02" value="103.9337005615234" />
    <PARAM id="df_room_decay_03" value="2.5" />
    <PARAM id="df_room_diffuse_03" value="69.20000457763672" />
    <PARAM id="df_room_early_03" value="-42.79999923706055" />
    <PARAM id="df_room_early_damp_03" value="98.96864318847656" />
    <PARAM id="df_room_early_send_03" value="-54.20000076293945" />
    <PARAM id="df_room_high_cut_03" value="123.3410873413086" />
    <PARAM id="df_room_late_03" value="-4.399997711181641" />
    <PARAM id="df_room_late_damp_03" value="104.9322357177734" />
    <PARAM id="df_room_low_boost_03" value="55.79999923706055" />
    <PARAM id="df_room_low_boost_freq_03" value="70.26498413085938" />
    <PARAM id="df_room_low_cut_03" value="-59.69720077514648" />
    <PARAM id="df_room_size_03" value="25.60000038146973" />
    <PARAM id="df_room_spin_03" value="-66.24396514892578" />
    <PARAM id="df_room_wander_03" value="28.39999961853027" />
    <PARAM id="echo_cycles_delay_beats_04" value="0.75" />
    <PARAM id="echo_cycles_feedback_ratio_04" value="57.19999694824219" />
    <PARAM id="echo_cycles_feedback_rotation_04" value="94.40000152587891" />
    <PARAM id="echo_cycles_filter_bandwidth_04" value="2.75" />
    <PARAM id="echo_cycles_filter_db_04" value="20.70000076293945" />
    <PARAM id="echo_cycles_filter_freq_04" value="82.16960906982422" />
    <PARAM id="echo_cycles_input_rotation_initial_04" value="88.40000152587891" />
    <PARAM id="echo_cycles_input_width_04" value="50.0" />
    <PARAM id="echo_cycles_output_variation_04" value="18.39999961853027" />
    <PARAM id="echo_cycles_rotation_mode_04" value="0.0" />
    <PARAM id="fx_mix_01" value="23.20000076293945" />
    <PARAM id="fx_mix_02" value="52.80000305175781" />
    <PARAM id="fx_mix_03" value="14.40000057220459" />
    <PARAM id="fx_mix_04" value="10.0" />
    <PARAM id="fx_type_01" value="27.0" />
    <PARAM id="fx_type_02" value="17.0" />
    <PARAM id="fx_type_03" value="33.0" />
    <PARAM id="fx_type_04" value="4.0" />
    <PARAM id="global_volume" value="-1.280000686645508" />
    <PARAM id="mixer_sends" value="7.0" />
    <PARAM id="wet_balance_02" value="0.0" />
    <PARAM id="wet_balance_03" value="-60.79999923706055" />
    <PARAM id="wet_pan_01" value="-55.20000076293945" />
    <PARAM id="wet_pan_02" value="55.19999694824219" />
    <PARAM id="wet_pan_03" value="44.80000305175781" />
    <PARAM id="wet_pan_04" value="-38.40000152587891" />
    <PARAM id="zebigchorus3_algo_01" value="4.0" />
    <PARAM id="zebigchorus3_delay_01" value="17.24844932556152" />
    <PARAM id="zebigchorus3_dispmod_01" value="0.4039999842643738" />
    <PARAM id="zebigchorus3_dispstatic_01" value="0.8279999494552612" />
    <PARAM id="zebigchorus3_dmod_01" value="68.0" />
    <PARAM id="zebigchorus3_g0_01" value="-0.1848000288009644" />
    <PARAM id="zebigchorus3_g1_01" value="-0.4444000124931335" />
    <PARAM id="zebigchorus3_grate_01" value="0.6459999680519104" />
    <PARAM id="zebigchorus3_gratedisp_01" value="0.6399999856948853" />
    <PARAM id="zebigchorus3_rate_01" value="0.1659999936819077" />
    <PARAM id="zebigchorus3_ratedisp_01" value="0.3240000009536743" />
    </params>
     )END"},
  {.name = "FX: Darkness",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="dry_pan_06" value="0.0" />
    <PARAM id="fdnverb_cascade_01" value="40.40000152587891" />
    <PARAM id="fdnverb_cascade_02" value="46.79999923706055" />
    <PARAM id="fdnverb_density_01" value="35.60000228881836" />
    <PARAM id="fdnverb_density_02" value="28.0" />
    <PARAM id="fdnverb_feedback_01" value="42.40000152587891" />
    <PARAM id="fdnverb_feedback_02" value="48.79999923706055" />
    <PARAM id="fdnverb_mod_depth_01" value="70.20000457763672" />
    <PARAM id="fdnverb_mod_depth_02" value="56.60000228881836" />
    <PARAM id="fdnverb_mod_rate_01" value="1.254000067710876" />
    <PARAM id="fdnverb_mod_rate_02" value="3.254000186920166" />
    <PARAM id="fdnverb_post_shift_01" value="-7.0" />
    <PARAM id="fdnverb_post_shift_02" value="-12.0" />
    <PARAM id="fdnverb_pre_shift_01" value="12.0" />
    <PARAM id="fdnverb_pre_shift_02" value="7.0" />
    <PARAM id="fdnverb_stereo_01" value="-29.40000152587891" />
    <PARAM id="fdnverb_stereo_02" value="-13.40000152587891" />
    <PARAM id="fdnverb_time_01" value="294.3999938964844" />
    <PARAM id="fdnverb_time_02" value="396.2000122070312" />
    <PARAM id="fx_mix_01" value="98.00000762939453" />
    <PARAM id="fx_mix_02" value="100.0" />
    <PARAM id="fx_mix_05" value="19.60000038146973" />
    <PARAM id="fx_mix_06" value="20.0" />
    <PARAM id="fx_mix_07" value="100.0" />
    <PARAM id="fx_type_01" value="38.0" />
    <PARAM id="fx_type_02" value="38.0" />
    <PARAM id="fx_type_03" value="0.0" />
    <PARAM id="fx_type_05" value="24.0" />
    <PARAM id="fx_type_06" value="37.0" />
    <PARAM id="global_volume" value="0.0" />
    <PARAM id="in_selection_02" value="1.0" />
    <PARAM id="in_selection_03" value="0.0" />
    <PARAM id="in_selection_04" value="1.0" />
    <PARAM id="in_selection_05" value="1.0" />
    <PARAM id="in_selection_06" value="0.0" />
    <PARAM id="in_selection_07" value="0.0" />
    <PARAM id="in_selection_08" value="2.0" />
    <PARAM id="mixer_sends" value="16.0" />
    <PARAM id="mute_solo_06" value="0.0" />
    <PARAM id="mute_solo_07" value="0.0" />
    <PARAM id="mute_solo_08" value="0.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="out_selection_03" value="0.0" />
    <PARAM id="out_selection_04" value="2.0" />
    <PARAM id="out_selection_05" value="0.0" />
    <PARAM id="out_selection_06" value="1.0" />
    <PARAM id="out_selection_07" value="0.0" />
    <PARAM id="out_selection_08" value="1.0" />
    <PARAM id="pan_07" value="0.0" />
    <PARAM id="routing" value="1.0" />
    <PARAM id="volume_05" value="-1.799999237060547" />
    <PARAM id="volume_06" value="-20.19000053405762" />
    <PARAM id="wet_pan_01" value="-35.20000457763672" />
    <PARAM id="wet_pan_02" value="38.39999389648438" />
    <PARAM id="wet_pan_05" value="0.0" />
    <PARAM id="wet_pan_07" value="0.0" />
    <PARAM id="chow_phaser_d1_05" value="0.1000000014901161" />
    <PARAM id="chow_phaser_d2_05" value="0.1000000014901161" />
    <PARAM id="chow_phaser_d3_05" value="0.1000000014901161" />
    <PARAM id="chow_phaser_feedback_05" value="0.8700000047683716" />
    <PARAM id="chow_phaser_freq_mult_05" value="0.0" />
    <PARAM id="chow_phaser_lfo_depth_05" value="0.7199999690055847" />
    <PARAM id="chow_phaser_lfo_freq_05" value="-14.76266479492188" />
    <PARAM id="chow_phaser_modulation_05" value="0.6999999284744263" />
    <PARAM id="chow_phaser_skew_05" value="1.420000076293945" />
    <PARAM id="chow_phaser_src_channel_05" value="0.0" />
    <PARAM id="chow_phaser_stages_05" value="40.91999816894531" />
    <PARAM id="fairly_childish_bias_05" value="56.0" />
    <PARAM id="fairly_childish_threshold_05" value="-31.69999885559082" />
    <PARAM id="fairly_childish_time_constant_05" value="2.0" />
    <PARAM id="filter2x_band1_drive_06" value="-5.799999237060547" />
    <PARAM id="filter2x_band1_feedback_06" value="71.19999694824219" />
    <PARAM id="filter2x_band1_freq_06" value="53.21419143676758" />
    <PARAM id="filter2x_band1_gain_06" value="-14.69999980926514" />
    <PARAM id="filter2x_band1_reso_06" value="0.6240000128746033" />
    <PARAM id="filter2x_band1_tolerance_06" value="-27.19999694824219" />
    <PARAM id="filter2x_band1_topology_06" value="1.0" />
    <PARAM id="filter2x_band1_type_06" value="1.0" />
    <PARAM id="filter2x_band2_drive_06" value="-30.39999961853027" />
    <PARAM id="filter2x_band2_feedback_06" value="68.80000305175781" />
    <PARAM id="filter2x_band2_freq_06" value="91.96116638183594" />
    <PARAM id="filter2x_band2_gain_06" value="2.700000762939453" />
    <PARAM id="filter2x_band2_reso_06" value="0.7719999551773071" />
    <PARAM id="filter2x_band2_tolerance_06" value="47.19999694824219" />
    <PARAM id="filter2x_band2_type_06" value="1.0" />
</params>
     )END"},
  {.name = "FX: Bloated Phaser.",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="atlantis_damping_03" value="33.0" />
    <PARAM id="atlantis_damping_high_03" value="123.9420928955078" />
    <PARAM id="atlantis_damping_low_03" value="63.15462112426758" />
    <PARAM id="atlantis_decay_03" value="1.490000009536743" />
    <PARAM id="atlantis_detune_shift_03" value="0.4160000085830688" />
    <PARAM id="atlantis_detune_spd_03" value="212.0" />
    <PARAM id="atlantis_shimmer_03" value="0.6299999952316284" />
    <PARAM id="atlantis_shimmer_tone_03" value="0.0" />
    <PARAM id="atlantis_window_03" value="40.59999847412109" />
    <PARAM id="bbd_echo_age_06" value="-55.20000076293945" />
    <PARAM id="bbd_echo_clock_curve_06" value="0.9399999976158142" />
    <PARAM id="bbd_echo_clock_offset_06" value="4.401000022888184" />
    <PARAM id="bbd_echo_clock_scale_06" value="0.4550000131130219" />
    <PARAM id="bbd_echo_delay_sync_l_06" value="28.0" />
    <PARAM id="bbd_echo_feedback_06" value="26.40000152587891" />
    <PARAM id="bbd_echo_hp_filter_06" value="29.0" />
    <PARAM id="bbd_echo_hp_res_06" value="0.8499999642372131" />
    <PARAM id="bbd_echo_lfo_depth_06" value="0.5600000023841858" />
    <PARAM id="bbd_echo_lfo_speed_06" value="0.3499999940395355" />
    <PARAM id="bbd_echo_lp_filter_06" value="95.0" />
    <PARAM id="bbd_echo_lp_res_06" value="0.949999988079071" />
    <PARAM id="bbd_echo_stages_06" value="3.0" />
    <PARAM id="channel_modifs_01" value="1.0" />
    <PARAM id="channel_modifs_02" value="8.0" />
    <PARAM id="eq4x_band1_freq_05" value="95.78802490234375" />
    <PARAM id="eq4x_band1_q_05" value="6.810199737548828" />
    <PARAM id="eq4x_band1_type_05" value="5.0" />
    <PARAM id="eq4x_band2_freq_05" value="55.30887985229492" />
    <PARAM id="eq4x_band2_q_05" value="9.626199722290039" />
    <PARAM id="eq4x_band2_type_05" value="6.0" />
    <PARAM id="eq4x_band3_freq_05" value="69.0" />
    <PARAM id="eq4x_band3_q_05" value="0.4999999701976776" />
    <PARAM id="eq4x_band3_type_05" value="0.0" />
    <PARAM id="fx_mix_01" value="62.40000152587891" />
    <PARAM id="fx_mix_02" value="34.40000152587891" />
    <PARAM id="fx_mix_03" value="100.0" />
    <PARAM id="fx_mix_05" value="100.0" />
    <PARAM id="fx_mix_06" value="32.40000152587891" />
    <PARAM id="fx_mix_07" value="11.60000038146973" />
    <PARAM id="fx_type_01" value="27.0" />
    <PARAM id="fx_type_02" value="3.0" />
    <PARAM id="fx_type_03" value="29.0" />
    <PARAM id="fx_type_05" value="23.0" />
    <PARAM id="fx_type_06" value="2.0" />
    <PARAM id="fx_type_07" value="25.0" />
    <PARAM id="global_volume" />
    <PARAM id="in_selection_02" value="0.0" />
    <PARAM id="in_selection_03" value="0.0" />
    <PARAM id="in_selection_04" value="1.0" />
    <PARAM id="in_selection_05" value="1.0" />
    <PARAM id="in_selection_06" value="0.0" />
    <PARAM id="in_selection_07" value="0.0" />
    <PARAM id="in_selection_08" value="2.0" />
    <PARAM id="mixer_sends" value="51.0" />
    <PARAM id="mute_solo_01" value="0.0" />
    <PARAM id="mute_solo_08" value="0.0" />
    <PARAM id="myphaser_feedback_07" value="66.40000915527344" />
    <PARAM id="myphaser_high_freq_07" value="80.4708251953125" />
    <PARAM id="myphaser_lfo_depth_07" value="31.0" />
    <PARAM id="myphaser_lfo_rate_07" value="-80.0" />
    <PARAM id="myphaser_lfo_rate_sync_07" value="24.0" />
    <PARAM id="myphaser_lfo_stereo_07" value="260.6000061035156" />
    <PARAM id="myphaser_lfo_wave_07" value="6.0" />
    <PARAM id="myphaser_low_freq_07" value="61.96545028686523" />
    <PARAM id="myphaser_q_07" value="0.1060000061988831" />
    <PARAM id="myphaser_stages_07" value="12.0" />
    <PARAM id="myphaser_stages_mode_07" value="6.0" />
    <PARAM id="out_selection_03" value="1.0" />
    <PARAM id="out_selection_04" value="2.0" />
    <PARAM id="out_selection_07" value="1.0" />
    <PARAM id="out_selection_08" value="1.0" />
    <PARAM id="routing" value="1.0" />
    <PARAM id="sandwitch_amp_secondary_gain_db_01" value="0.0" />
    <PARAM id="sandwitch_amp_secondary_gain_db_03" value="-25.59999847412109" />
    <PARAM id="sound_delay_l_02" value="0.0" />
    <PARAM id="sound_delay_l_03" value="0.0" />
    <PARAM id="sound_delay_l_sync_02" value="1.0" />
    <PARAM id="sound_delay_r_02" value="0.0" />
    <PARAM id="sound_delay_r_03" value="0.0" />
    <PARAM id="sound_delay_r_sync_02" value="2.0" />
    <PARAM id="volume_01" value="0.0" />
    <PARAM id="volume_03" value="-3.960002899169922" />
    <PARAM id="volume_05" value="6.599998474121094" />
    <PARAM id="volume_06" value="-5.180000305175781" />
    <PARAM id="volume_07" value="-5.430000305175781" />
    <PARAM id="volume_08" value="0.0" />
    <PARAM id="wet_pan_01" value="-11.20000457763672" />
    <PARAM id="wet_pan_02" value="25.59999847412109" />
    <PARAM id="wet_pan_03" value="-41.60000228881836" />
    <PARAM id="wet_pan_05" value="33.59999084472656" />
    <PARAM id="wet_pan_06" value="-28.0" />
    <PARAM id="wet_pan_07" value="24.79999542236328" />
    <PARAM id="zebigchorus3_algo_01" value="2.0" />
    <PARAM id="zebigchorus3_delay_01" value="17.0" />
    <PARAM id="zebigchorus3_dispstatic_01" value="0.4479999840259552" />
    <PARAM id="zebigchorus3_dmod_01" value="22.79999923706055" />
    <PARAM id="zebigchorus3_g1_01" value="-0.1556000113487244" />
    <PARAM id="zebigchorus3_grate_01" value="0.3059999942779541" />
    <PARAM id="zebigchorus3_gratedisp_01" value="0.699999988079071" />
    <PARAM id="zebigchorus3_rate_01" value="0.1579999923706055" />
    <PARAM id="zebigchorus3_ratedisp_01" value="0.3599999845027924" />
</params>
    )END"},
  {.name = "FX: Multiband Dynamics",
   .xml = R"END(
?xml version="1.0" encoding="UTF-8"?><params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="in_selection_02" value="1.0" />
    <PARAM id="in_selection_03" value="1.0" />
    <PARAM id="in_selection_04" value="1.0" />
    <PARAM id="in_selection_05" value="1.0" />
    <PARAM id="in_selection_06" value="2.0" />
    <PARAM id="in_selection_07" value="4.0" />
    <PARAM id="in_selection_08" value="8.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="out_selection_02" value="2.0" />
    <PARAM id="out_selection_03" value="4.0" />
    <PARAM id="out_selection_04" value="8.0" />
    <PARAM id="out_selection_05" value="1.0" />
    <PARAM id="out_selection_06" value="1.0" />
    <PARAM id="out_selection_07" value="1.0" />
    <PARAM id="out_selection_08" value="1.0" />
    <PARAM id="channel_modifs_01" value="0.0" />
    <PARAM id="channel_modifs_02" value="0.0" />
    <PARAM id="eq4x_band1_freq_01" value="57.99777221679688" />
    <PARAM id="eq4x_band1_freq_02" value="57.99777221679688" />
    <PARAM id="eq4x_band1_freq_03" value="84.78579711914062" />
    <PARAM id="eq4x_band1_freq_04" value="108.2253265380859" />
    <PARAM id="eq4x_band1_freq_05" value="69.0" />
    <PARAM id="eq4x_band1_q_01" value="0.4999999701976776" />
    <PARAM id="eq4x_band1_q_05" value="0.4999999701976776" />
    <PARAM id="eq4x_band1_type_01" value="5.0" />
    <PARAM id="eq4x_band1_type_02" value="6.0" />
    <PARAM id="eq4x_band1_type_03" value="6.0" />
    <PARAM id="eq4x_band1_type_04" value="6.0" />
    <PARAM id="eq4x_band1_type_05" value="0.0" />
    <PARAM id="eq4x_band2_freq_02" value="84.78579711914062" />
    <PARAM id="eq4x_band2_freq_03" value="108.2253265380859" />
    <PARAM id="eq4x_band2_freq_05" value="69.0" />
    <PARAM id="eq4x_band2_q_05" value="0.4999999701976776" />
    <PARAM id="eq4x_band2_type_02" value="5.0" />
    <PARAM id="eq4x_band2_type_03" value="5.0" />
    <PARAM id="eq4x_band2_type_05" value="0.0" />
    <PARAM id="fx_mix_01" value="100.0" />
    <PARAM id="fx_mix_02" value="100.0" />
    <PARAM id="fx_mix_06" value="100.0" />
    <PARAM id="fx_mix_07" value="78.80000305175781" />
    <PARAM id="fx_mix_08" value="74.80000305175781" />
    <PARAM id="fx_type_01" value="23.0" />
    <PARAM id="fx_type_02" value="23.0" />
    <PARAM id="fx_type_03" value="23.0" />
    <PARAM id="fx_type_04" value="23.0" />
    <PARAM id="fx_type_05" value="15.0" />
    <PARAM id="fx_type_06" value="39.0" />
    <PARAM id="fx_type_07" value="39.0" />
    <PARAM id="fx_type_08" value="11.0" />
    <PARAM id="gate_expander_gatethresh_08" value="-0.5" />
    <PARAM id="gate_expander_hysteresis_08" value="-3.100000381469727" />
    <PARAM id="gate_expander_range_08" value="-33.91999816894531" />
    <PARAM id="gate_expander_sidechain_freq_08" value="15.48682022094727" />
    <PARAM id="global_volume" value="-2.450000762939453" />
    <PARAM id="mixer_sends" value="0.0" />
    <PARAM id="mute_solo_01" value="0.0" />
    <PARAM id="mute_solo_02" value="0.0" />
    <PARAM id="mute_solo_03" value="0.0" />
    <PARAM id="mute_solo_05" value="0.0" />
    <PARAM id="mute_solo_06" value="0.0" />
    <PARAM id="mute_solo_07" value="0.0" />
    <PARAM id="mute_solo_08" value="0.0" />
    <PARAM id="routing" value="1.0" />
    <PARAM id="track_comp_autogain_06" value="1.0" />
    <PARAM id="track_comp_autogain_07" value="1.0" />
    <PARAM id="track_comp_compattack_06" value="3.560999870300293" />
    <PARAM id="track_comp_compattack_07" value="50.95100021362305" />
    <PARAM id="track_comp_compfeedbk_06" value="0.0" />
    <PARAM id="track_comp_compfeedbk_07" value="40.59999847412109" />
    <PARAM id="track_comp_compknee_07" value="5.170000076293945" />
    <PARAM id="track_comp_comprange_07" value="-28.95999908447266" />
    <PARAM id="track_comp_compratio_07" value="7.700000286102295" />
    <PARAM id="track_comp_comprelease_06" value="235.25" />
    <PARAM id="track_comp_comprelease_07" value="74.8699951171875" />
    <PARAM id="track_comp_compthresh_06" value="-20.15999984741211" />
    <PARAM id="track_comp_compthresh_07" value="-39.59999847412109" />
    <PARAM id="track_comp_compwindow_07" value="9.590000152587891" />
    <PARAM id="track_comp_linkamount_07" value="9.0" />
    <PARAM id="track_comp_saturation_06" value="0.0" />
    <PARAM id="track_comp_saturation_07" value="0.0" />
    <PARAM id="track_comp_scfreq_06" value="20.0" />
    <PARAM id="track_comp_scfreq_07" value="20.0" />
    <PARAM id="transience_attack_05" value="0.0" />
    <PARAM id="transience_attack_amt_05" value="-1.0" />
    <PARAM id="transience_decay_05" value="0.3919999897480011" />
    <PARAM id="transience_decay_amt_05" value="0.7839999198913574" />
    <PARAM id="transience_gainsmoothing_05" value="0.1799999922513962" />
    <PARAM id="transience_mode_05" value="0.0" />
    <PARAM id="volume_03" value="0.0" />
    <PARAM id="volume_05" value="0.3400001525878906" />
    <PARAM id="volume_06" value="3.149997711181641" />
    <PARAM id="volume_07" value="2.219997406005859" />
    <PARAM id="volume_08" value="5.919998168945312" />
    <PARAM id="wet_pan_07" value="-18.40000152587891" />
    <PARAM id="wet_pan_08" value="48.0" />
</params>
     )END"},
  {.name = "FX: Multiband Nasal Modulations",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="in_selection_02" value="0.0" />
    <PARAM id="in_selection_03" value="0.0" />
    <PARAM id="in_selection_04" value="0.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="out_selection_03" value="1.0" />
    <PARAM id="out_selection_04" value="1.0" />
    <PARAM id="fx_type_01" value="42.0" />
    <PARAM id="fx_type_02" value="27.0" />
    <PARAM id="fx_type_03" value="2.0" />
    <PARAM id="fx_type_04" value="25.0" />
    <PARAM id="bbd_echo_age_07" value="0.0" />
    <PARAM id="bbd_echo_clock_scale_07" value="1.0" />
    <PARAM id="bbd_echo_delay_sync_l_07" value="50.0" />
    <PARAM id="bbd_echo_delay_sync_r_07" value="50.0" />
    <PARAM id="bbd_echo_feedback_07" value="40.0" />
    <PARAM id="bbd_echo_stages_07" value="9.0" />
    <PARAM id="fx_mix_06" value="100.0" />
    <PARAM id="fx_mix_07" value="100.0" />
    <PARAM id="fx_mix_08" value="100.0" />
    <PARAM id="fx_type_05" value="0.0" />
    <PARAM id="fx_type_06" value="0.0" />
    <PARAM id="fx_type_07" value="0.0" />
    <PARAM id="fx_type_08" value="0.0" />
    <PARAM id="global_volume" value="0.0" />
    <PARAM id="in_selection_05" value="0.0" />
    <PARAM id="in_selection_06" value="0.0" />
    <PARAM id="in_selection_07" value="0.0" />
    <PARAM id="in_selection_08" value="0.0" />
    <PARAM id="mixer_sends" value="4.0" />
    <PARAM id="routing" value="0.0" />
    <PARAM id="bbd_echo_age_03" value="-74.40000152587891" />
    <PARAM id="bbd_echo_clock_scale_03" value="0.8490000367164612" />
    <PARAM id="bbd_echo_delay_sync_l_03" value="25.59999847412109" />
    <PARAM id="bbd_echo_delay_sync_r_03" value="18.39999961853027" />
    <PARAM id="bbd_echo_feedback_03" value="56.0" />
    <PARAM id="bbd_echo_stages_03" value="2.0" />
    <PARAM id="lr_crossv_band1_diff" value="-11.20000457763672" />
    <PARAM id="lr_crossv_band1_frequency" value="57.78551864624023" />
    <PARAM id="lr_crossv_band1_mode" />
    <PARAM id="lr_crossv_band1_out" value="0.0" />
    <PARAM id="lr_crossv_band2_diff" value="24.0" />
    <PARAM id="lr_crossv_band2_frequency" value="84.61087036132812" />
    <PARAM id="lr_crossv_band2_mode" />
    <PARAM id="lr_crossv_band2_out" value="1.0" />
    <PARAM id="lr_crossv_band3_diff" value="-19.20000457763672" />
    <PARAM id="lr_crossv_band3_frequency" value="106.9496917724609" />
    <PARAM id="lr_crossv_band3_mode" />
    <PARAM id="lr_crossv_band3_out" value="2.0" />
    <PARAM id="fx_mix_02" value="76.80000305175781" />
    <PARAM id="fx_mix_03" value="57.60000228881836" />
    <PARAM id="fx_mix_04" value="53.60000228881836" />
    <PARAM id="myphaser_feedback_04" value="-74.40000152587891" />
    <PARAM id="myphaser_high_freq_04" value="107.3255462646484" />
    <PARAM id="myphaser_lfo_depth_04" value="65.80000305175781" />
    <PARAM id="myphaser_lfo_rate_04" value="-80.0" />
    <PARAM id="myphaser_lfo_rate_sync_04" value="16.0" />
    <PARAM id="myphaser_lfo_stereo_04" value="298.1000061035156" />
    <PARAM id="myphaser_low_freq_04" value="84.93194580078125" />
    <PARAM id="myphaser_q_04" value="1.632000088691711" />
    <PARAM id="myphaser_stages_04" value="11.0" />
    <PARAM id="myphaser_stages_mode_04" value="4.0" />
    <PARAM id="volume_03" value="2.919998168945312" />
    <PARAM id="volume_04" value="2.219997406005859" />
    <PARAM id="wet_balance_02" value="5.860000610351562" />
    <PARAM id="wet_balance_06" value="0.0" />
    <PARAM id="wet_pan_02" value="54.39999389648438" />
    <PARAM id="wet_pan_03" value="-34.40000152587891" />
    <PARAM id="wet_pan_04" value="-59.20000076293945" />
    <PARAM id="zebigchorus3_algo_02" value="7.0" />
    <PARAM id="zebigchorus3_delay_02" value="17.79999923706055" />
</params>
     )END"},
  {.name = "FX: Earthquake",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="bbe_cv_04" value="8.800000190734863" />
    <PARAM id="bbe_drive_04" value="11.0" />
    <PARAM id="bbe_low_contour_04" value="0.0" />
    <PARAM id="bbe_process_04" value="26.80000114440918" />
    <PARAM id="fx_mix_02" value="100.0" />
    <PARAM id="fx_mix_03" value="60.00000381469727" />
    <PARAM id="fx_mix_04" value="24.40000152587891" />
    <PARAM id="fx_mix_05" value="100.0" />
    <PARAM id="fx_mix_06" value="100.0" />
    <PARAM id="fx_mix_07" value="14.80000114440918" />
    <PARAM id="fx_mix_08" value="6.800000190734863" />
    <PARAM id="fx_type_01" value="42.0" />
    <PARAM id="fx_type_02" value="11.0" />
    <PARAM id="fx_type_03" value="15.0" />
    <PARAM id="fx_type_04" value="20.0" />
    <PARAM id="fx_type_05" value="41.0" />
    <PARAM id="fx_type_06" value="41.0" />
    <PARAM id="fx_type_07" value="27.0" />
    <PARAM id="fx_type_08" value="34.0" />
    <PARAM id="gate_expander_attack_02" value="7.700000286102295" />
    <PARAM id="gate_expander_gatethresh_02" value="-17.0" />
    <PARAM id="gate_expander_hysteresis_02" value="-4.440000057220459" />
    <PARAM id="gate_expander_range_02" value="-24.15999984741211" />
    <PARAM id="gate_expander_release_02" value="678.5" />
    <PARAM id="global_volume" value="-3.649999618530273" />
    <PARAM id="in_selection_05" value="1.0" />
    <PARAM id="in_selection_06" value="2.0" />
    <PARAM id="in_selection_07" value="4.0" />
    <PARAM id="in_selection_08" value="8.0" />
    <PARAM id="lr_crossv_band1_diff" value="0.0" />
    <PARAM id="lr_crossv_band1_frequency" />
    <PARAM id="lr_crossv_band1_mode" value="3.0" />
    <PARAM id="lr_crossv_band1_out" value="0.0" />
    <PARAM id="lr_crossv_band2_diff" value="-2.400001525878906" />
    <PARAM id="lr_crossv_band2_frequency" value="77.43550872802734" />
    <PARAM id="lr_crossv_band2_mode" />
    <PARAM id="lr_crossv_band2_out" value="1.0" />
    <PARAM id="lr_crossv_band3_diff" value="3.199996948242188" />
    <PARAM id="lr_crossv_band3_frequency" value="106.9496917724609" />
    <PARAM id="lr_crossv_band3_mode" />
    <PARAM id="lr_crossv_band3_out" value="2.0" />
    <PARAM id="mixer_sends" value="0.0" />
    <PARAM id="mute_solo_01" value="0.0" />
    <PARAM id="mute_solo_02" value="0.0" />
    <PARAM id="mute_solo_03" value="0.0" />
    <PARAM id="mute_solo_05" value="0.0" />
    <PARAM id="mute_solo_06" value="0.0" />
    <PARAM id="mute_solo_07" value="0.0" />
    <PARAM id="mute_solo_08" value="0.0" />
    <PARAM id="out_selection_02" value="2.0" />
    <PARAM id="out_selection_03" value="4.0" />
    <PARAM id="out_selection_04" value="8.0" />
    <PARAM id="out_selection_05" value="1.0" />
    <PARAM id="out_selection_06" value="1.0" />
    <PARAM id="out_selection_07" value="1.0" />
    <PARAM id="out_selection_08" value="1.0" />
    <PARAM id="routing" value="1.0" />
    <PARAM id="tal_reverb2_decay_08" value="13.60000038146973" />
    <PARAM id="tal_reverb2_lowshelf_frequency_08" value="114.347412109375" />
    <PARAM id="tal_reverb2_stereo_width_08" value="83.20000457763672" />
    <PARAM id="volume_01" value="0.0" />
    <PARAM id="volume_03" value="4.30999755859375" />
    <PARAM id="volume_04" value="1.520000457763672" />
    <PARAM id="volume_06" value="0.0" />
    <PARAM id="volume_07" value="0.0" />
    <PARAM id="volume_08" value="-1.799999237060547" />
    <PARAM id="waveshaper_band_gain_05" value="-2.75" />
    <PARAM id="waveshaper_band_gain_06" value="0.0" />
    <PARAM id="waveshaper_drive_05" value="-4.5" />
    <PARAM id="waveshaper_drive_06" value="-16.25" />
    <PARAM id="waveshaper_drive_balance_06" value="-13.75" />
    <PARAM id="waveshaper_emphasis_amount_05" value="1.25" />
    <PARAM id="waveshaper_emphasis_amount_06" value="20.75" />
    <PARAM id="waveshaper_emphasis_freq_05" value="38.36567687988281" />
    <PARAM id="waveshaper_emphasis_freq_06" value="73.51918029785156" />
    <PARAM id="waveshaper_emphasis_q_05" value="0.4999999701976776" />
    <PARAM id="waveshaper_emphasis_q_06" value="5.160000324249268" />
    <PARAM id="waveshaper_envfollow_attack_05" value="14.30000019073486" />
    <PARAM id="waveshaper_envfollow_attack_06" value="19.30000114440918" />
    <PARAM id="waveshaper_envfollow_mode_05" value="3.0" />
    <PARAM id="waveshaper_envfollow_mode_06" value="1.0" />
    <PARAM id="waveshaper_envfollow_release_05" value="480.0" />
    <PARAM id="waveshaper_envfollow_release_06" value="481.0" />
    <PARAM id="waveshaper_envfollow_sensitivity_05" value="-10.23999977111816" />
    <PARAM id="waveshaper_envfollow_sensitivity_06" value="20.0" />
    <PARAM id="waveshaper_envfollow_to_dc_06" value="0.0" />
    <PARAM id="waveshaper_envfollow_to_drive_05" value="0.0" />
    <PARAM id="waveshaper_envfollow_to_drive_06" value="100.0" />
    <PARAM id="waveshaper_envfollow_to_emphasis_amount_05" value="0.0" />
    <PARAM id="waveshaper_envfollow_to_emphasis_amount_06" value="-92.80000305175781" />
    <PARAM id="waveshaper_envfollow_to_emphasis_freq_05" value="-60.0" />
    <PARAM id="waveshaper_envfollow_to_emphasis_freq_06" value="-90.40000152587891" />
    <PARAM id="waveshaper_feedback_05" value="86.10000610351562" />
    <PARAM id="waveshaper_feedback_06" value="4.599998474121094" />
    <PARAM id="waveshaper_mode_06" value="2.0" />
    <PARAM id="waveshaper_type_05" value="1.0" />
    <PARAM id="waveshaper_type_06" value="2.0" />
    <PARAM id="wet_balance_01" value="0.0" />
    <PARAM id="wet_balance_02" value="0.0" />
    <PARAM id="wet_balance_04" value="0.0" />
    <PARAM id="wet_balance_08" value="0.0" />
    <PARAM id="wet_pan_03" value="0.0" />
    <PARAM id="wet_pan_07" value="-44.0" />
    <PARAM id="wet_pan_08" value="48.80000305175781" />
    <PARAM id="zebigchorus3_algo_07" value="5.0" />
    <PARAM id="zebigchorus3_delay_07" value="7.799999713897705" />
    <PARAM id="zebigchorus3_dispmod_07" value="0.8919999599456787" />
    <PARAM id="zebigchorus3_dispstatic_07" value="0.8399999737739563" />
    <PARAM id="zebigchorus3_dmod_07" value="92.39999389648438" />
    <PARAM id="zebigchorus3_g0_07" value="0.4459999203681946" />
    <PARAM id="zebigchorus3_g1_07" value="-0.4824000000953674" />
    <PARAM id="zebigchorus3_grate_07" value="0.1099999994039536" />
    <PARAM id="zebigchorus3_gratedisp_07" value="0.7519999742507935" />
    <PARAM id="zebigchorus3_rate_07" value="0.1180000007152557" />
    <PARAM id="zebigchorus3_ratedisp_07" value="0.1319999992847443" />
</params>
     )END"},
  {.name = "FX: Chorused Reverb",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1005000">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="fx_type_01" value="0.0" />
    <PARAM id="fx_type_02" value="50.0" />
    <PARAM id="in_selection_02" value="1.0" />
    <PARAM id="in_selection_03" value="0.0" />
    <PARAM id="mixer_sends" value="0.0" />
    <PARAM id="out_selection_02" value="1.0" />
    <PARAM id="out_selection_03" value="0.0" />
    <PARAM id="reverb_damp_factor_02" value="10.39999961853027" />
    <PARAM id="reverb_damp_freq_02" value="-33.59999847412109" />
    <PARAM id="reverb_decay_02" value="13401.0" />
    <PARAM id="reverb_diff_in_02" value="0.0" />
    <PARAM id="reverb_diff_out_02" value="22.0" />
    <PARAM id="reverb_early_gain_02" value="-0.3999977111816406" />
    <PARAM id="reverb_early_size_02" value="-75.19999694824219" />
    <PARAM id="reverb_gap_02" value="3.000000238418579" />
    <PARAM id="reverb_hp_freq_02" value="35.0" />
    <PARAM id="reverb_l_sparseness_02" value="-31.26000213623047" />
    <PARAM id="reverb_late_gain_02" value="-6.0" />
    <PARAM id="reverb_lf_time_factor_02" value="0.0" />
    <PARAM id="reverb_lr_sparseness_02" value="11.62999725341797" />
    <PARAM id="reverb_mod_depth_02" value="89.19999694824219" />
    <PARAM id="reverb_mod_freq_02" value="47.20000076293945" />
    <PARAM id="reverb_mod_mode_02" value="1.0" />
    <PARAM id="reverb_mod_spread_02" value="100.0" />
    <PARAM id="reverb_predelay_02" value="1.0" />
    <PARAM id="reverb_r_sparseness_02" value="77.89999389648438" />
    <PARAM id="reverb_size_02" value="43.20001220703125" />
    <PARAM id="reverb_stereo_02" value="100.0" />
    <PARAM id="track_comp_compattack_01" value="10.00100040435791" />
    <PARAM id="track_comp_compfeedbk_01" value="25.0" />
    <PARAM id="track_comp_compratio_01" value="4.0" />
    <PARAM id="track_comp_comprelease_01" value="200.0" />
    <PARAM id="track_comp_compthresh_01" value="0.0" />
    <PARAM id="track_comp_saturation_01" value="25.0" />
    <PARAM id="volume_02" value="-1.08" />
    <PARAM id="wet_balance_03" value="-5.860000610351562" />
    <PARAM id="wet_pan_02" value="14.39999389648438" />
    <gui_data bus1_text="Dry" bus2_text="Wet" bus3_text="" bus4_text="" bus5_text="" bus6_text="" bus7_text="" bus8_text="" />
</params>
     )END"},
  {.name = "FX: Shimmer Reverb",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1005000">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="fx_type_01" value="49.0" />
    <PARAM id="fx_type_02" value="49.0" />
    <PARAM id="in_selection_02" value="1.0" />
    <PARAM id="in_selection_03" value="1.0" />
    <PARAM id="mixer_sends" value="16.0" />
    <PARAM id="out_selection_02" value="2.0" />
    <PARAM id="out_selection_03" value="4.0" />
    <PARAM id="volume_02" value="-6.900001525878906" />
    <PARAM id="wet_balance_03" value="-5.860000610351562" />
    <PARAM id="wet_pan_02" value="14.39999389648438" />
    <PARAM id="filter2x_band1_ef_to_freq_05" value="79.19999694824219" />
    <PARAM id="filter2x_band1_freq_05" value="76.65371704101562" />
    <PARAM id="filter2x_band1_reso_05" value="0.7519999742507935" />
    <PARAM id="filter2x_band1_type_05" value="1.0" />
    <PARAM id="filter2x_envfollow_attack_05" value="27.20000076293945" />
    <PARAM id="filter2x_envfollow_release_05" value="227.0" />
    <PARAM id="filter2x_envfollow_sensitivity_05" value="27.04000091552734" />
    <PARAM id="fx_mix_02" value="100.0" />
    <PARAM id="fx_type_03" value="0.0" />
    <PARAM id="fx_type_05" value="37.0" />
    <PARAM id="fx_type_06" value="50.0" />
    <PARAM id="global_volume" value="0.4399986267089844" />
    <PARAM id="in_selection_05" value="3.0" />
    <PARAM id="in_selection_07" value="4.0" />
    <PARAM id="naive_pitch_amt_06" value="12.0" />
    <PARAM id="out_selection_06" value="1.0" />
    <PARAM id="out_selection_07" value="1.0" />
    <PARAM id="reverb_damp_factor_06" value="7.599999904632568" />
    <PARAM id="reverb_damp_freq_06" value="0.0" />
    <PARAM id="reverb_decay_06" value="17152.0" />
    <PARAM id="reverb_diff_out_06" value="12.79999923706055" />
    <PARAM id="reverb_early_gain_06" value="-0.3999977111816406" />
    <PARAM id="reverb_early_size_06" value="-91.19999694824219" />
    <PARAM id="reverb_l_sparseness_06" value="-31.26000213623047" />
    <PARAM id="reverb_lf_time_factor_06" value="-4.800003051757812" />
    <PARAM id="reverb_lr_sparseness_06" value="11.62999725341797" />
    <PARAM id="reverb_mod_depth_06" value="68.0" />
    <PARAM id="reverb_mod_freq_06" value="55.87000274658203" />
    <PARAM id="reverb_mod_mode_06" value="3.0" />
    <PARAM id="reverb_mod_spread_06" value="40.0" />
    <PARAM id="reverb_r_sparseness_06" value="77.89999389648438" />
    <PARAM id="reverb_size_06" value="-100.0" />
    <PARAM id="routing" value="1.0" />
    <PARAM id="rubberband_formant_mode_01" value="0.0" />
    <PARAM id="rubberband_mode_01" value="0.0" />
    <PARAM id="rubberband_semitones_01" value="12.0" />
    <PARAM id="soundtouch_detune_02" value="-0.1200000047683716" />
    <PARAM id="soundtouch_mode_01" value="2.0" />
    <PARAM id="soundtouch_mode_02" value="1.0" />
    <PARAM id="soundtouch_semitones_01" value="-12.0" />
    <PARAM id="soundtouch_semitones_02" value="12.0" />
    <PARAM id="volume_01" value="0.0" />
    <PARAM id="volume_05" value="0.80999755859375" />
    <PARAM id="volume_06" value="-4.94000244140625" />
    <PARAM id="volume_07" value="0.0" />
    <PARAM id="wet_pan_03" value="-5.860000610351562" />
    <gui_data bus1_text="-12" bus2_text="+12" bus3_text="Dry Passthrough" bus4_text="" bus5_text="Pre-filter" bus6_text="Wet" bus7_text="Dry" bus8_text="" />
</params>
     )END"},
  {.name = "FX: Phaser DDL-ish (Drum/full mix)",
   .xml = R"END(
<?xml version="1.0" encoding="UTF-8"?>
<params plugin_version="1009000">
    <PARAM id="in_selection_01" value="1.0" />
    <PARAM id="out_selection_01" value="1.0" />
    <PARAM id="fx_type_01" value="43.0" />
    <PARAM id="fx_type_02" value="25.0" />
    <PARAM id="in_selection_02" value="0.0" />
    <PARAM id="in_selection_03" value="0.0" />
    <PARAM id="mixer_sends" value="2.0" />
    <PARAM id="out_selection_02" value="0.0" />
    <PARAM id="out_selection_03" value="1.0" />
    <PARAM id="diffuse_delay_damp_03" value="0.0" />
    <PARAM id="diffuse_delay_desync_03" value="64.79999542236328" />
    <PARAM id="diffuse_delay_diffusion_03" value="0.0" />
    <PARAM id="diffuse_delay_feedback_03" value="32.53000259399414" />
    <PARAM id="diffuse_delay_freq_spread_03" value="84.0" />
    <PARAM id="diffuse_delay_hipass_03" value="43.0" />
    <PARAM id="diffuse_delay_mod_depth_03" value="1.209999918937683" />
    <PARAM id="diffuse_delay_mod_freq_03" value="0.1410000026226044" />
    <PARAM id="diffuse_delay_mod_mode_03" value="2.0" />
    <PARAM id="diffuse_delay_mode_03" value="2.0" />
    <PARAM id="diffuse_delay_sixteenths_03" value="8.029999732971191" />
    <PARAM id="diffuse_delay_transients_03" value="-100.0" />
    <PARAM id="dry_balance_02" value="0.0" />
    <PARAM id="dry_balance_03" value="19.19999694824219" />
    <PARAM id="fx_mix_03" value="8.0" />
    <PARAM id="fx_type_03" value="51.0" />
    <PARAM id="global_volume" value="0.0" />
    <PARAM id="in_selection_04" value="1.0" />
    <PARAM id="mute_solo_03" value="0.0" />
    <PARAM id="myphaser_delay_lfo_02" value="24.00000762939453" />
    <PARAM id="myphaser_delay_time_02" value="5.600000381469727" />
    <PARAM id="myphaser_feedback_02" value="100.0" />
    <PARAM id="myphaser_high_freq_02" value="123.0762329101562" />
    <PARAM id="myphaser_lfo_depth_02" value="100.0" />
    <PARAM id="myphaser_lfo_rate_02" value="-70.22529602050781" />
    <PARAM id="myphaser_lfo_stereo_02" value="31.60000038146973" />
    <PARAM id="myphaser_lfo_wave_02" value="1.0" />
    <PARAM id="myphaser_low_freq_02" value="74.9683837890625" />
    <PARAM id="myphaser_parallel_mix_02" value="-74.40000152587891" />
    <PARAM id="myphaser_q_02" value="1.304000020027161" />
    <PARAM id="myphaser_stages_02" value="10.0" />
    <PARAM id="myphaser_stages_mode_02" value="12.0" />
    <PARAM id="myphaser_topology_02" value="3.0" />
    <PARAM id="out_selection_04" value="1.0" />
    <PARAM id="pan_01" value="0.0" />
    <PARAM id="volume_01" value="0.0" />
    <PARAM id="volume_02" value="0.0" />
    <PARAM id="volume_03" value="2.919998168945312" />
    <PARAM id="volume_04" value="-1.5" />
    <PARAM id="volume_07" value="0.0" />
    <PARAM id="wet_balance_01" value="0.0" />
    <PARAM id="wonky_crossv_band1_frequency" value="41.31813049316406" />
    <PARAM id="wonky_crossv_band1_mode" value="1.0" />
    <PARAM id="wonky_crossv_band1_out" value="0.0" />
    <PARAM id="wonky_crossv_band2_mode" value="0.0" />
    <PARAM id="wonky_crossv_band2_out" value="0.0" />
    <PARAM id="wonky_crossv_band3_mode" value="0.0" />
    <PARAM id="wonky_crossv_band3_out" value="0.0" />
    <gui_data bus1_text=" Wet Lows" bus2_text="Pass -&gt;" bus3_text="Wet Highs" bus4_text="Dry" bus5_text="" bus6_text="" bus7_text="" bus8_text="" />
</params>
     )END"},
};

// clang-format on
#if 0
  {.name = "Routing: ",
   .xml = R"END(

     )END"},
#endif
// REMINDER. unused parameter cleanup regex: .*_[\d][\d]" />\n
} // namespace artv
