#include "mix-maxtrix/order_and_buffering.hpp"
#include "gtest/gtest.h"

namespace artv {

TEST (order_and_buffering, ins_and_outs)
{
  order_and_buffering<4> ob;
  ob.recompute (0x8421, 0x8421, 0, 0, 0, {});

  ASSERT_EQ (ob.in[0][0], 1);
  ASSERT_EQ (ob.in[0][1], 0);
  ASSERT_EQ (ob.in[0][2], 0);
  ASSERT_EQ (ob.in[0][3], 0);

  ASSERT_EQ (ob.in[1][0], 0);
  ASSERT_EQ (ob.in[1][1], 2);
  ASSERT_EQ (ob.in[1][2], 0);
  ASSERT_EQ (ob.in[1][3], 0);

  ASSERT_EQ (ob.in[2][0], 0);
  ASSERT_EQ (ob.in[2][1], 0);
  ASSERT_EQ (ob.in[2][2], 3);
  ASSERT_EQ (ob.in[2][3], 0);

  ASSERT_EQ (ob.in[3][0], 0);
  ASSERT_EQ (ob.in[3][1], 0);
  ASSERT_EQ (ob.in[3][2], 0);
  ASSERT_EQ (ob.in[3][3], 4);

  ASSERT_EQ (ob.out[0][0], 1);
  ASSERT_EQ (ob.out[0][1], 0);
  ASSERT_EQ (ob.out[0][2], 0);
  ASSERT_EQ (ob.out[0][3], 0);

  ASSERT_EQ (ob.out[1][0], 0);
  ASSERT_EQ (ob.out[1][1], 2);
  ASSERT_EQ (ob.out[1][2], 0);
  ASSERT_EQ (ob.out[1][3], 0);

  ASSERT_EQ (ob.out[2][0], 0);
  ASSERT_EQ (ob.out[2][1], 0);
  ASSERT_EQ (ob.out[2][2], 3);
  ASSERT_EQ (ob.out[2][3], 0);

  ASSERT_EQ (ob.out[3][0], 0);
  ASSERT_EQ (ob.out[3][1], 0);
  ASSERT_EQ (ob.out[3][2], 0);
  ASSERT_EQ (ob.out[3][3], 4);
}

TEST (order_and_buffering, passthrough_order)
{
  order_and_buffering<4> ob;
  ob.recompute (0x8421, 0x8421, 0, 0, 0, {});
  ASSERT_EQ (ob.order[0], 0);
  ASSERT_EQ (ob.order[1], 1);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);
}

TEST (order_and_buffering, input_dependency)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x9421, // channel 1 enabled as input on channel 4
    0x8421,
    0,
    0,
    0,
    {});
  ASSERT_EQ (ob.order[0], 1);
  ASSERT_EQ (ob.order[1], 2);
  ASSERT_EQ (ob.order[2], 3);
  ASSERT_EQ (ob.order[3], 0); // channel 1 is processed last
}

TEST (order_and_buffering, input_dependency_no_modification)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x8431, // channel 1 enabled as input on channel 4
    0x8421,
    0,
    0,
    0,
    {},
    0,
    order_and_buffering<4>::bus_latency_arr {},
    [] (uint i) { return false; });
  ASSERT_EQ (ob.order[0], 0);
  ASSERT_EQ (ob.order[1], 1);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);
}

TEST (order_and_buffering, output_dependency)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x8421,
    0x8439, // own mix to out + mix2 on out1 +  mix1 on out 4
    0,
    0,
    0,
    {});
  ASSERT_EQ (ob.order[0], 1);
  ASSERT_EQ (ob.order[1], 2);
  ASSERT_EQ (ob.order[2], 3);
  ASSERT_EQ (ob.order[3], 0); // chnl1 after chnl4, as its out will sum.
}

TEST (order_and_buffering, output_dependency_no_modification)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x8421,
    0x8429, // own mix to out +  mix1 on out 4
    0,
    0,
    0,
    {});
  // chnl 1 kept on its position, it doesn't destroy its out buffer
  ASSERT_EQ (ob.order[0], 0);
  ASSERT_EQ (ob.order[1], 1);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);
}

TEST (order_and_buffering, mixed_dependencies)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x8423, // channel 1 has inputs from channel 1 and 2.
    0x8439, // own mix to out + mix2 on out1 +  mix1 on out 4
    0,
    0,
    0,
    {});
  ASSERT_EQ (ob.order[0], 2);
  ASSERT_EQ (ob.order[1], 3);
  ASSERT_EQ (ob.order[2], 0);
  ASSERT_EQ (ob.order[3], 1);
}

TEST (order_and_buffering, bad_case)
{
  order_and_buffering<4> ob;
  ob.recompute (
    // in1   int2     in3     in4
    0x0111 | 0x0222 | 0x0444 | 0x8000,
    // out1  out2     out3     out4
    0x0111 | 0x0222 | 0x0444 | 0x8008,
    0,
    0,
    0,
    {});
  ASSERT_EQ (ob.order[0], 3);
  ASSERT_EQ (ob.order[1], 0);
  ASSERT_EQ (ob.order[2], 1);
  ASSERT_EQ (ob.order[3], 2);
}

TEST (order_and_buffering, out_modification_by_no_output)
{
  order_and_buffering<4> ob;
  ob.recompute (0x8421, 0x8422, 0, 0, 0, {});
  ASSERT_EQ (ob.order[0], 1);
  ASSERT_EQ (ob.order[1], 0);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);
}

TEST (order_and_buffering, in_modification_by_no_input)
{
  order_and_buffering<4> ob;
  ob.recompute (0x8430, 0x8421, 0, 0, 0, {});
  ASSERT_EQ (ob.order[0], 1);
  ASSERT_EQ (ob.order[1], 0);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);
}

TEST (order_and_buffering, pre_channel_processing_swap)
{
  order_and_buffering<4> ob;
  ob.recompute (0x1428, 0x8421, 0, 0, 0, {});
  ASSERT_EQ (ob.swaps[0], 4);
  ASSERT_EQ (ob.swaps[1], 0);
  ASSERT_EQ (ob.swaps[2], 0);
  ASSERT_EQ (ob.swaps[3], 0);

  ASSERT_EQ (ob.order[0], 0);
  ASSERT_EQ (ob.order[1], 1);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);

  ASSERT_EQ (ob.in[0][0], 1);
  ASSERT_EQ (ob.in[3][3], 4);
}

TEST (order_and_buffering, buffering_with_no_deps)
{
  order_and_buffering<4> ob;
  ob.recompute (0x8421, 0x8421, 0, 0, 0, {});
  ASSERT_EQ (ob.order[0], 0);
  ASSERT_EQ (ob.order[1], 1);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);

  ASSERT_EQ (ob.swaps[0], 0);
  ASSERT_EQ (ob.swaps[1], 0);
  ASSERT_EQ (ob.swaps[2], 0);
  ASSERT_EQ (ob.swaps[3], 0);

  ASSERT_EQ (ob.mix[0], 1);
  ASSERT_EQ (ob.mix[1], 2);
  ASSERT_EQ (ob.mix[2], 3);
  ASSERT_EQ (ob.mix[3], 4);
}

TEST (order_and_buffering, in_buffering_with_in_deps)
{
  order_and_buffering<4> ob;
  ob.recompute (0xb433, 0x8421, 0, 0, 0, {});

  ASSERT_EQ (ob.order[0], 2);
  ASSERT_EQ (ob.order[1], 3);
  ASSERT_EQ (ob.order[2], 0);
  ASSERT_EQ (ob.order[3], 1);

  ASSERT_EQ (ob.swaps[0], 0);
  ASSERT_EQ (ob.swaps[1], 0);
  ASSERT_EQ (ob.swaps[2], 0);
  ASSERT_EQ (ob.swaps[3], 0);

  ASSERT_EQ (ob.mix[0], -1);
  ASSERT_EQ (ob.mix[1], 2);
  ASSERT_EQ (ob.mix[2], 3);
  ASSERT_EQ (ob.mix[3], 4);
}

TEST (order_and_buffering, in_buffering_with_in_deps2)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x7777, // all depend on channels 1, 2 and 3
    0x8421,
    0,
    0,
    0,
    {});

  ASSERT_EQ (ob.order[0], 3);
  ASSERT_EQ (ob.order[1], 0);
  ASSERT_EQ (ob.order[2], 1);
  ASSERT_EQ (ob.order[3], 2);

  ASSERT_EQ (ob.swaps[0], 0);
  ASSERT_EQ (ob.swaps[1], 0);
  ASSERT_EQ (ob.swaps[2], 0);
  ASSERT_EQ (ob.swaps[3], 0);

  ASSERT_EQ (ob.mix[0], -1);
  ASSERT_EQ (ob.mix[1], -2);
  ASSERT_EQ (ob.mix[2], 3);
  ASSERT_EQ (ob.mix[3], 4);
}

TEST (order_and_buffering, out_buffering_with_in_deps)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x8421,
    // out1  out2     out3     out4
    0x0011 | 0x0022 | 0x0400 | 0x8088,
    0,
    0,
    0,
    {});

  ASSERT_EQ (ob.order[0], 2);
  ASSERT_EQ (ob.order[1], 3);
  ASSERT_EQ (ob.order[2], 0);
  ASSERT_EQ (ob.order[3], 1);

  ASSERT_EQ (ob.swaps[0], 0);
  ASSERT_EQ (ob.swaps[1], 0);
  ASSERT_EQ (ob.swaps[2], 0);
  ASSERT_EQ (ob.swaps[3], 0);

  ASSERT_EQ (ob.mix[0], -1);
  ASSERT_EQ (ob.mix[1], 2);
  ASSERT_EQ (ob.mix[2], 3);
  ASSERT_EQ (ob.mix[3], 4);
}

TEST (order_and_buffering, out_buffering_with_in_dep2)
{
  order_and_buffering<4> ob;
  ob.recompute (
    0x8421,
    // all depend on channels 1, 2 and 3
    // out1  out2     out3     out4
    0x0111 | 0x0222 | 0x0444 | 0x0888,
    0,
    0,
    0,
    {});

  ASSERT_EQ (ob.order[0], 3);
  ASSERT_EQ (ob.order[1], 0);
  ASSERT_EQ (ob.order[2], 1);
  ASSERT_EQ (ob.order[3], 2);

  ASSERT_EQ (ob.swaps[0], 0);
  ASSERT_EQ (ob.swaps[1], 0);
  ASSERT_EQ (ob.swaps[2], 0);
  ASSERT_EQ (ob.swaps[3], 0);

  ASSERT_EQ (ob.mix[0], -1);
  ASSERT_EQ (ob.mix[1], -2);
  ASSERT_EQ (ob.mix[2], 3);
  ASSERT_EQ (ob.mix[3], 4);
}

TEST (order_and_buffering, solo)
{
  order_and_buffering<4> ob;
  ob.recompute (0x8421, 0x8421, 0x0002, 0, 0, {});
  ASSERT_EQ (ob.in[0][0], 1);
  ASSERT_EQ (ob.in[1][1], 0);
  ASSERT_EQ (ob.in[2][2], 0);
  ASSERT_EQ (ob.in[3][3], 0);

  ASSERT_EQ (ob.out[0][0], 1);
  ASSERT_EQ (ob.out[1][1], 0);
  ASSERT_EQ (ob.out[2][2], 0);
  ASSERT_EQ (ob.out[3][3], 0);
}

TEST (order_and_buffering, mute)
{
  order_and_buffering<4> ob;
  constexpr auto         mute_chnl = 4;
  ob.recompute (0x8421, 0x8421, 1 << ((mute_chnl - 1) * 2), 0, 0, {});
  ASSERT_EQ (ob.in[0][0], 1);
  ASSERT_EQ (ob.in[1][1], 2);
  ASSERT_EQ (ob.in[2][2], 3);
  ASSERT_EQ (ob.in[3][3], 0);

  ASSERT_EQ (ob.out[0][0], 1);
  ASSERT_EQ (ob.out[1][1], 2);
  ASSERT_EQ (ob.out[2][2], 3);
  ASSERT_EQ (ob.out[3][3], 0);
}

TEST (order_and_buffering, mute_8)
{
  order_and_buffering<8> ob;
  constexpr auto         mute_chnl = 6;
  ob.recompute (
    0x8040201008040201,
    0x8040201008040201,
    1 << ((mute_chnl - 1) * 2),
    0,
    0,
    {});
  ASSERT_EQ (ob.in[0][0], 1);
  ASSERT_EQ (ob.in[1][1], 2);
  ASSERT_EQ (ob.in[2][2], 3);
  ASSERT_EQ (ob.in[3][3], 4);
  ASSERT_EQ (ob.in[4][4], 5);
  ASSERT_EQ (ob.in[5][5], 0);
  ASSERT_EQ (ob.in[6][6], 7);
  ASSERT_EQ (ob.in[7][7], 8);

  ASSERT_EQ (ob.out[0][0], 1);
  ASSERT_EQ (ob.out[1][1], 2);
  ASSERT_EQ (ob.out[2][2], 3);
  ASSERT_EQ (ob.out[3][3], 4);
  ASSERT_EQ (ob.out[4][4], 5);
  ASSERT_EQ (ob.out[5][5], 0);
  ASSERT_EQ (ob.out[6][6], 7);
  ASSERT_EQ (ob.out[7][7], 8);
}

TEST (order_and_buffering, mixer2mixer_sends_inhibit_reordering)
{
  order_and_buffering<4> ob;
  ob.recompute (0x8431, 0x8421, 0, 1, 0, {});
  // As in2 depends on in1, normally mix1 would be processed after mix2, but
  // as mix1 is sent to mix2, the reordering is inhibited for mix1 and 2.
  ASSERT_EQ (ob.order[0], 0);
  ASSERT_EQ (ob.order[1], 1);
  ASSERT_EQ (ob.order[2], 2);
  ASSERT_EQ (ob.order[3], 3);

  ASSERT_EQ (ob.swaps[0], 0);
  ASSERT_EQ (ob.swaps[1], 0);
  ASSERT_EQ (ob.swaps[2], 0);
  ASSERT_EQ (ob.swaps[3], 0);

  ASSERT_EQ (ob.mix[0], -1);
  ASSERT_EQ (ob.mix[1], 2);
  ASSERT_EQ (ob.mix[2], 3);
  ASSERT_EQ (ob.mix[3], 4);

  ASSERT_EQ (ob.in[1][0], 1);
  ASSERT_EQ (ob.in[1][1], 2);
  ASSERT_EQ (ob.in[1][2], 0);
  ASSERT_EQ (ob.in[1][3], 0);

  ASSERT_EQ (ob.receives[1][0], -1);

  // regular stuff
  ASSERT_EQ (ob.in[0][0], 1);
  ASSERT_EQ (ob.in[1][1], 2);
  ASSERT_EQ (ob.in[2][2], 3);
  ASSERT_EQ (ob.in[3][3], 4);

  ASSERT_EQ (ob.out[0][0], -1);
  ASSERT_EQ (ob.out[1][1], 2);
  ASSERT_EQ (ob.out[2][2], 3);
  ASSERT_EQ (ob.out[3][3], 4);
}

TEST (order_and_buffering, latency_empty)
{
  order_and_buffering<4>::bus_latency_arr lat = {0, 0, 0, 0};
  order_and_buffering<4>                  ob;
  ob.recompute (0x1111, 0x1111, 0, 0, 0, {}, 0, lat);

  ASSERT_EQ (ob.receives_latency[0], 0);
  ASSERT_EQ (ob.receives_latency[1], 0);
  ASSERT_EQ (ob.receives_latency[2], 0);
  ASSERT_EQ (ob.receives_latency[3], 0);

  ASSERT_EQ (ob.pre_output_mix_latency[0], 0);
  ASSERT_EQ (ob.pre_output_mix_latency[1], 0);
  ASSERT_EQ (ob.pre_output_mix_latency[2], 0);
  ASSERT_EQ (ob.pre_output_mix_latency[3], 0);

  ASSERT_EQ (ob.plugin_latency, 0);
}

TEST (order_and_buffering, latency_no_sends)
{
  order_and_buffering<4>::bus_latency_arr lat = {1, 2, 3, 4};
  order_and_buffering<4>                  ob;
  ob.recompute (0x1111, 0x1111, 0, 0, 0, {}, 0, lat);

  ASSERT_EQ (ob.receives_latency[0], 0);
  ASSERT_EQ (ob.receives_latency[1], 0);
  ASSERT_EQ (ob.receives_latency[2], 0);
  ASSERT_EQ (ob.receives_latency[3], 0);

  ASSERT_EQ (ob.pre_output_mix_latency[0], 3);
  ASSERT_EQ (ob.pre_output_mix_latency[1], 2);
  ASSERT_EQ (ob.pre_output_mix_latency[2], 1);
  ASSERT_EQ (ob.pre_output_mix_latency[3], 0);

  ASSERT_EQ (ob.plugin_latency, 4);
}

TEST (order_and_buffering, latency_all_sends_enabled)
{
  order_and_buffering<4>::bus_latency_arr lat = {1, 2, 3, 4};
  order_and_buffering<4>                  ob;
  ob.recompute (0x1111, 0x1111, 0, 0x7, 0, {}, 0, lat);

  ASSERT_EQ (ob.receives_latency[0], 0);
  ASSERT_EQ (ob.receives_latency[1], 1);
  ASSERT_EQ (ob.receives_latency[2], 3);
  ASSERT_EQ (ob.receives_latency[3], 6);

  // cummulative latencies are [1, 3, 6, 10]

  ASSERT_EQ (ob.pre_output_mix_latency[0], 10 - 1);
  ASSERT_EQ (ob.pre_output_mix_latency[1], 10 - 3);
  ASSERT_EQ (ob.pre_output_mix_latency[2], 10 - 6);
  ASSERT_EQ (ob.pre_output_mix_latency[3], 10 - 10);

  ASSERT_EQ (ob.plugin_latency, 10);
}

TEST (order_and_buffering, latency_all_sends_enabled_series)
{
  order_and_buffering<4>::bus_latency_arr lat = {1, 2, 3, 4};
  order_and_buffering<4>                  ob;
  // buses 1 to 4 connected in series. Only bus 1 inputs and bus 4 outputs
  ob.recompute (0x0001, 0x1000, 0, 0x7, 0, {}, 0, lat);

  ASSERT_EQ (ob.fx_latency[0], 1);
  ASSERT_EQ (ob.fx_latency[1], 2);
  ASSERT_EQ (ob.fx_latency[2], 3);
  ASSERT_EQ (ob.fx_latency[3], 4);

  ASSERT_EQ (ob.receives_latency[0], 0);
  ASSERT_EQ (ob.receives_latency[1], 0);
  ASSERT_EQ (ob.receives_latency[2], 0);
  ASSERT_EQ (ob.receives_latency[3], 0);

  // cummulative latencies are [1, 3, 6, 10]

  ASSERT_EQ (ob.pre_output_mix_latency[0], 0);
  ASSERT_EQ (ob.pre_output_mix_latency[1], 0);
  ASSERT_EQ (ob.pre_output_mix_latency[2], 0);
  ASSERT_EQ (ob.pre_output_mix_latency[3], 0);

  ASSERT_EQ (ob.plugin_latency, 10);
}

TEST (order_and_buffering, latency_no_sends_groups_enabled)
{
  order_and_buffering<4>::bus_latency_arr lat = {1, 2, 3, 4};
  order_and_buffering<4>                  ob;
  ob.recompute (0x1111, 0x1111, 0, 0, 0, {}, 1, lat);

  ASSERT_EQ (ob.fx_latency[0], 1);
  ASSERT_EQ (ob.fx_latency[1], 2);
  ASSERT_EQ (ob.fx_latency[2], 3);
  ASSERT_EQ (ob.fx_latency[3], 4);

  ASSERT_EQ (ob.receives_latency[0], 0);
  ASSERT_EQ (ob.receives_latency[1], 0);
  ASSERT_EQ (ob.receives_latency[2], 0);
  ASSERT_EQ (ob.receives_latency[3], 0);

  ASSERT_EQ (ob.pre_output_mix_latency[0], 2 - 1);
  ASSERT_EQ (ob.pre_output_mix_latency[1], 2 - 2);
  ASSERT_EQ (ob.pre_output_mix_latency[2], 4 - 3);
  ASSERT_EQ (ob.pre_output_mix_latency[3], 4 - 4);

  ASSERT_EQ (ob.plugin_latency, 2 + 4);
}

TEST (order_and_buffering, latency_all_sends_enabled_groups_enabled)
{
  order_and_buffering<4>::bus_latency_arr lat = {1, 2, 3, 4};
  order_and_buffering<4>                  ob;
  ob.recompute (0x1111, 0x1111, 0, 0x7, 0, {}, 1, lat);

  ASSERT_EQ (ob.fx_latency[0], 1);
  ASSERT_EQ (ob.fx_latency[1], 2);
  ASSERT_EQ (ob.fx_latency[2], 3);
  ASSERT_EQ (ob.fx_latency[3], 4);

  ASSERT_EQ (ob.receives_latency[0], 0);
  ASSERT_EQ (ob.receives_latency[1], 1);
  ASSERT_EQ (ob.receives_latency[2], 0);
  ASSERT_EQ (ob.receives_latency[3], 3);

  // cummulative latencies are [1, 3], [3, 7]

  ASSERT_EQ (ob.pre_output_mix_latency[0], 3 - 1);
  ASSERT_EQ (ob.pre_output_mix_latency[1], 3 - 3);
  ASSERT_EQ (ob.pre_output_mix_latency[2], 7 - 3);
  ASSERT_EQ (ob.pre_output_mix_latency[3], 7 - 7);

  ASSERT_EQ (ob.plugin_latency, 3 + 7);
}

} // namespace artv
