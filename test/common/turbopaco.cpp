#include <array>
#include <gtest/gtest.h>
#include <vector>

#include "artv-common/dsp/own/fx/turbopaco-engine.hpp"

namespace artv {

static constexpr uint max_block_size_test = 8;
static constexpr uint delay_size
  = max_block_size_test * 2 - max_block_size_test / 2;

struct delay_test_spec {
  static constexpr auto get_spec()
  {
    return make_array<detail::turbopaco::stage_data> (
      detail::turbopaco::make_block_delay (delay_size, 1));
  }
};

//------------------------------------------------------------------------------
TEST (turbopaco, block_processed_delays)
{
  using namespace detail::turbopaco;

  engine<delay_test_spec, delay::data_type::float32, max_block_size_test> e;

  std::vector<u8>                        mem;
  std::array<float, max_block_size_test> in, out;

  mem.resize (e.get_required_bytes());
  e.reset_memory (mem);

  s16 in_counter  = 0;
  s16 out_counter = 0;
  // block 1
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 0);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  for (uint i = 0; i < max_block_size_test; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // block 2
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 0);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  auto expected_data_bound = max_block_size_test / 2;
  for (uint i = 0; i < expected_data_bound; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // first sample should come out
  for (uint i = expected_data_bound; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 3
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 0);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  for (uint i = 0; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 4
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 0);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  for (uint i = 0; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
}
//------------------------------------------------------------------------------
TEST (turbopaco, block_processed_delays_witch_shift)
{
  using namespace detail::turbopaco;

  engine<delay_test_spec, delay::data_type::float32, max_block_size_test> e;
  std::vector<u8>                                                         mem;
  std::array<float, max_block_size_test> in, out;

  mem.resize (e.get_required_bytes());
  e.reset_memory (mem);

  s16 in_counter  = 0;
  s16 out_counter = 0;
  // block 1
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 1);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  for (uint i = 0; i < max_block_size_test; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // block 2
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 1);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  auto expected_data_bound = (max_block_size_test / 2) + 1;
  for (uint i = 0; i < expected_data_bound; ++i) {
    EXPECT_EQ (0.f, out[i]);
  }
  // first sample should come out
  for (uint i = expected_data_bound; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 3
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 1);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  for (uint i = 0; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
  // block 4
  for (uint i = 0; i < max_block_size_test; ++i) {
    in[i] = ((float) ++in_counter) * 0.01;
  }
  e.fetch_block (stage_list<0> {}, xspan {out}, 1);
  e.push (stage_list<0> {}, xspan {in}.to_const());

  for (uint i = 0; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, out[i], 0.001f);
  }
}
//------------------------------------------------------------------------------
TEST (turbopaco, plain_delays)
{
  using namespace detail::turbopaco;
  engine<delay_test_spec, delay::data_type::float32, max_block_size_test> e;
  std::vector<u8>                                                         mem;
  std::array<float, max_block_size_test>                                  io;

  mem.resize (e.get_required_bytes());
  e.reset_memory (mem);

  s16 in_counter  = 0;
  s16 out_counter = 0;
  // block 1
  for (uint i = 0; i < max_block_size_test; ++i) {
    io[i] = ((float) ++in_counter) * 0.01;
  }
  e.run (stage_list<0> {}, xspan {io});
  for (uint i = 0; i < max_block_size_test; ++i) {
    EXPECT_EQ (0.f, io[i]);
  }
  // block 2
  for (uint i = 0; i < max_block_size_test; ++i) {
    io[i] = ((float) ++in_counter) * 0.01;
  }
  e.run (stage_list<0> {}, xspan {io});
  auto expected_data_bound = (max_block_size_test / 2);
  for (uint i = 0; i < expected_data_bound; ++i) {
    EXPECT_EQ (0.f, io[i]);
  }
  // first sample should come out
  for (uint i = expected_data_bound; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, io[i], 0.001f);
  }
  // block 3
  for (uint i = 0; i < max_block_size_test; ++i) {
    io[i] = ((float) ++in_counter) * 0.01;
  }
  e.run (stage_list<0> {}, xspan {io});
  for (uint i = 0; i < max_block_size_test; ++i) {
    ++out_counter;
    EXPECT_NEAR (((float) out_counter) * 0.01, io[i], 0.001f);
  }
}
//------------------------------------------------------------------------------
} // namespace artv
