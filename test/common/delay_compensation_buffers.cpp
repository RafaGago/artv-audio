#include <array>
#include <gtest/gtest.h>

#include "artv-common/misc/delay_compensation_buffers.hpp"

namespace artv {
//------------------------------------------------------------------------------
TEST (block_delay_compensation, bigger_input)
{
  block_delay_compensation<uint, 1> dly;

  using in_arr = std::array<uint, 4>;
  std::array<uint, 1> latencies {2};

  dly.reset<uint> (make_crange (latencies));

  in_arr in {1, 2, 3, 4};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 0);
  ASSERT_EQ (in[1], 0);
  ASSERT_EQ (in[2], 1);
  ASSERT_EQ (in[3], 2);

  in = in_arr {5, 6, 7, 8};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 3);
  ASSERT_EQ (in[1], 4);
  ASSERT_EQ (in[2], 5);
  ASSERT_EQ (in[3], 6);

  in = in_arr {9, 10, 11, 12};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 7);
  ASSERT_EQ (in[1], 8);
  ASSERT_EQ (in[2], 9);
  ASSERT_EQ (in[3], 10);
}
//------------------------------------------------------------------------------
TEST (block_delay_compensation, bigger_delay)
{
  block_delay_compensation<uint, 1> dly;

  using in_arr = std::array<uint, 2>;
  std::array<uint, 1> latencies {4};

  dly.reset<uint> (make_crange (latencies));

  in_arr in {1, 2};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 0);
  ASSERT_EQ (in[1], 0);

  in = in_arr {3, 4};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 0);
  ASSERT_EQ (in[1], 0);

  in = in_arr {5, 6};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 1);
  ASSERT_EQ (in[1], 2);

  in = in_arr {7, 8};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 3);
  ASSERT_EQ (in[1], 4);

  in = in_arr {9, 10};
  dly.compensate (0, in);

  ASSERT_EQ (in[0], 5);
  ASSERT_EQ (in[1], 6);
}
//------------------------------------------------------------------------------
} // namespace artv
