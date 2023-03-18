#include <array>
#include <gtest/gtest.h>
#include <vector>

#include "artv-common/dsp/own/fx/lofiverb-engine.hpp"

namespace artv {

static constexpr uint max_block_size = 8;
static constexpr uint delay_size     = max_block_size * 2 - max_block_size / 2;

static constexpr auto get_delay_test_spec()
{
  return make_array<detail::lofiverb::stage_data> (
    detail::lofiverb::make_block_delay (delay_size, 1));
}

struct delay_test_spec {
  static constexpr auto values {get_delay_test_spec()};
};

//------------------------------------------------------------------------------
TEST (lofiverb, block_processed_delays)
{
  detail::lofiverb::engine<delay_test_spec, max_block_size> e;
  std::vector<s16>                                          mem;
  std::array<float, max_block_size>                         in, out;

  mem.resize (e.get_required_size());
  e.reset_memory (mem);

  s16 in_counter  = 0;
  s16 out_counter = 0;
  // block 1
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out});
  e.push<0> (xspan {in}.to_const());

  for (uint i = 0; i < max_block_size; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // block 2
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out});
  e.push<0> (xspan {in}.to_const());

  auto expected_data_bound = max_block_size / 2;
  for (uint i = 0; i < expected_data_bound; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // first sample should come out
  for (uint i = expected_data_bound; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 3
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out});
  e.push<0> (xspan {in}.to_const());

  for (uint i = 0; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 4
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out});
  e.push<0> (xspan {in}.to_const());

  for (uint i = 0; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
}
//------------------------------------------------------------------------------
TEST (lofiverb, block_processed_delays_witch_shift)
{
  detail::lofiverb::engine<delay_test_spec, max_block_size> e;
  std::vector<s16>                                          mem;
  std::array<float, max_block_size>                         in, out;

  mem.resize (e.get_required_size());
  e.reset_memory (mem);

  s16 in_counter  = 0;
  s16 out_counter = 0;
  // block 1
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out}, 1);
  e.push<0> (xspan {in}.to_const());

  for (uint i = 0; i < max_block_size; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // block 2
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out}, 1);
  e.push<0> (xspan {in}.to_const());

  auto expected_data_bound = (max_block_size / 2) + 1;
  for (uint i = 0; i < expected_data_bound; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // first sample should come out
  for (uint i = expected_data_bound; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 3
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out}, 1);
  e.push<0> (xspan {in}.to_const());

  for (uint i = 0; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 4
  for (uint i = 0; i < max_block_size; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block<0> (xspan {out}, 1);
  e.push<0> (xspan {in}.to_const());

  for (uint i = 0; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
}
//------------------------------------------------------------------------------
TEST (lofiverb, plain_delays)
{
  detail::lofiverb::engine<delay_test_spec, max_block_size> e;
  std::vector<s16>                                          mem;
  std::array<float, max_block_size>                         io;

  mem.resize (e.get_required_size());
  e.reset_memory (mem);

  s16 in_counter  = 0;
  s16 out_counter = 0;
  // block 1
  for (uint i = 0; i < max_block_size; ++i) {
    io[i] = ((float) ++in_counter) * 0.01;
  }
  e.run<0> (xspan {io});
  for (uint i = 0; i < max_block_size; ++i) {
    EXPECT_EQ (0.f, io[i]);
  }
  // block 2
  for (uint i = 0; i < max_block_size; ++i) {
    io[i] = ((float) ++in_counter) * 0.01;
  }
  e.run<0> (xspan {io});
  auto expected_data_bound = (max_block_size / 2);
  for (uint i = 0; i < expected_data_bound; ++i) {
    EXPECT_EQ (0.f, io[i]);
  }
  // first sample should come out
  for (uint i = expected_data_bound; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, io[i], 0.001f);
  }
  // block 3
  for (uint i = 0; i < max_block_size; ++i) {
    io[i] = ((float) ++in_counter) * 0.01;
  }
  e.run<0> (xspan {io});
  for (uint i = 0; i < max_block_size; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, io[i], 0.001f);
  }
}
//------------------------------------------------------------------------------
} // namespace artv
