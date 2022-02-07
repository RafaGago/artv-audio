// Sanity testing FFTs to see if they match what jsfx expects

#include <cstring>

#include <math.h>

#include "gtest/gtest.h"

#include "artv-common/dsp/own/classes/mufft.hpp"
#include "artv-common/dsp/own/classes/wdl_fft.hpp"

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

//------------------------------------------------------------------------------
template <class T, class FFT>
class fft_test : public ::testing::Test {
public:
  fft_test() {}
  //----------------------------------------------------------------------------
  static constexpr uint block_size = 16;
  //----------------------------------------------------------------------------
  void SetUp()
  {
    in_mem.resize (block_size * 2 / simdwrapper::n_elems);
    out_mem.resize (block_size * 2 / simdwrapper::n_elems);

    in  = &in_mem[0][0];
    out = &out_mem[0][0];

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

  using simdwrapper = simd_mem<T, 32 / sizeof (T), 32>;

  std::vector<simdwrapper> in_mem {};
  std::vector<simdwrapper> out_mem {};

  T* in  = nullptr;
  T* out = nullptr;
};
//------------------------------------------------------------------------------
using mufft_test  = fft_test<float, artv::mufft::initialized_ffts<float, 16>>;
using wdlfft_test = fft_test<double, artv::wdl::initialized_ffts<double>>;
//------------------------------------------------------------------------------
TEST_F (wdlfft_test, wdl_matches_in)
{
  fft_impl.forward_transform (out, in, block_size);
  fft_impl.backward_transform (out, block_size);
  scale_buffer (out);
  match_fft_buffers (out, in);
}
//------------------------------------------------------------------------------
TEST_F (wdlfft_test, permuted_wdl_matches_in)
{
  fft_impl.forward_transform (out, in, block_size);
  fft_impl.forward_permute (out, block_size);
  fft_impl.backward_permute (out, block_size);
  fft_impl.backward_transform (out, block_size);
  scale_buffer (out);
  match_fft_buffers (out, in);
}
//------------------------------------------------------------------------------
TEST_F (mufft_test, mufft_matches_in)
{
  fft_impl.forward_transform (out, in, block_size);
  fft_impl.backward_transform (out, block_size);
  scale_buffer (out);
  match_fft_buffers (out, in);
}
//------------------------------------------------------------------------------
} // namespace artv
