// Sanity testing FFTs to see if they match what jsfx expects

#include <cstring>
#include <math.h>
#include <vector>

#include "gtest/gtest.h"

#include "artv-common/dsp/own/classes/fft.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"

#if 0
// JSFX verification
in = 0;
block_size = 16;

increment=(8*$pi) / (block_size - 1);
i=0;

loop(block_size,
    in[i * 2] = sin (i * increment);
    in[i * 2 + 1] = 0;
    i += 1;
    );

in0 = in[0];
in1 = in[2];
in2 = in[4];
in3 = in[6];
in4 = in[8];
in5 = in[10];
in6 = in[12];
in7 = in[14];
in8 = in[16];
in9 = in[18];
in10 = in[20];
in11 = in[22];
in12 = in[24];
in13 = in[26];
in14 = in[28];
in15 = in[30];

fft(in, block_size);
fft_permute(in, block_size);
fft_ipermute(in, block_size);
ifft(in, block_size);

i=0;
loop(block_size,
    in[i * 2] /= block_size;
    in[i * 2 + 1] /= block_size;
    i += 1;
    );

out0 = in[0];
out1 = in[2];
out2 = in[4];
out3 = in[6];
out4 = in[8];
out5 = in[10];
out6 = in[12];
out7 = in[14];
out8 = in[16];
out9 = in[18];
out10 = in[20];
out11 = in[22];
out12 = in[24];
out13 = in[26];
out14 = in[28];
out15 = in[30];

out_i_0 = in[1];
out_i_1 = in[3];
out_i_2 = in[5];
out_i_3 = in[7];
out_i_4 = in[9];
out_i_5 = in[11];
out_i_6 = in[13];
out_i_7 = in[15];
out_i_8 = in[17];
out_i_9 = in[19];
out_i_10 = in[21];
out_i_11 = in[23];
out_i_12 = in[25];
out_i_13 = in[27];
out_i_14 = in[29];
out_i_15 = in[31];

#endif

namespace artv {

static constexpr uint fftsize = 512;

//------------------------------------------------------------------------------
template <class T, class FFT>
class fft_test : public ::testing::Test {
public:
  fft_test() {}
  //----------------------------------------------------------------------------
  static constexpr uint block_size = fftsize;
  //----------------------------------------------------------------------------
  void SetUp()
  {
    in.resize (block_size * 2);
    out.resize (block_size * 2);

    constexpr double pi = M_PI;

    double increment = (8. * pi) / (block_size - 1);

    for (uint i = 0; i < block_size; ++i) {
      in[i * 2]     = (T) (sin (double (i) * increment));
      in[i * 2 + 1] = 0.f;
    }
  }
  //----------------------------------------------------------------------------
  void TearDown() {}
  //----------------------------------------------------------------------------
  static void match_fft_buffers (T* a, T* b)
  {
    for (uint i = 0; i < block_size; ++i) {
      ASSERT_NEAR (a[i * 2], b[i * 2], 0.000001);
      ASSERT_NEAR (a[i * 2 + 1], b[i * 2 + 1], 0.000001);
    }
  }
  //----------------------------------------------------------------------------
  void scale_buffer (T* b)
  {
    float scale = 1. / ((double) block_size);
    for (uint i = 0; i < block_size * 2; ++i) {
      b[i] *= scale;
    }
  }
  //----------------------------------------------------------------------------
  FFT fft_impl;
  std::vector<T, overaligned_allocator<T, FFT::fft_type::io_alignment>> in;
  std::vector<T, overaligned_allocator<T, FFT::fft_type::io_alignment>> out;
};
//------------------------------------------------------------------------------
using mufft_test
  = fft_test<float, artv::initialized_ffts<float, true, fftsize>>;
using pffft_test
  = fft_test<double, artv::initialized_ffts<double, true, fftsize>>;
//------------------------------------------------------------------------------
TEST_F (pffft_test, pffft_matches_in)
{
  auto& fft = *fft_impl.get_fft (block_size);
  fft.forward (out, in);
  fft.backward (out);
  scale_buffer (out.data());
  match_fft_buffers (out.data(), in.data());
}
//------------------------------------------------------------------------------
TEST_F (pffft_test, permuted_pffft_matches_in)
{
  auto& fft = *fft_impl.get_fft (block_size);
  fft.forward (out, in);
  fft.reorder_after_forward (out);
  fft.reorder_before_backward (out);
  fft.backward (out);
  scale_buffer (out.data());
  match_fft_buffers (out.data(), in.data());
}
//------------------------------------------------------------------------------
TEST_F (mufft_test, mufft_matches_in)
{
  auto& fft = *fft_impl.get_fft (block_size);
  fft.forward (out, in);
  fft.backward (out);
  scale_buffer (out.data());
  match_fft_buffers (out.data(), in.data());
}
//------------------------------------------------------------------------------
} // namespace artv
