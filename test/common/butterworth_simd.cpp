#include <array>
#include <cstdio>
#include <cstring>
#include <gtest/gtest.h>

#include "artv-common/dsp/own/blocks/filters/composite/butterworth.hpp"
#include "common/simd_block_test.hpp"

namespace artv {

//------------------------------------------------------------------------------
TEST (butterworth, lowpass_7_poles)
{
  using dsp      = butterworth<7>;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f = 500.;

  dsp::lowpass (test.scalar_coeffs, f, test.samplerate);
  dsp::lowpass_simd (
    test.vector_coeffs, vec_set<svf_test::vector_type> (f), test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (butterworth, highpass_7_poles)
{
  using dsp      = butterworth<7>;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f = 500.;

  dsp::highpass ({test.scalar_coeffs}, f, test.samplerate);
  dsp::highpass_simd (
    test.vector_coeffs, vec_set<svf_test::vector_type> (f), test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (butterworth, lowpass_8_poles)
{
  using dsp      = butterworth<8>;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f = 500.;

  dsp::lowpass ({test.scalar_coeffs}, f, test.samplerate);
  dsp::lowpass_simd (
    test.vector_coeffs, vec_set<svf_test::vector_type> (f), test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (butterworth, highpass_8_poles)
{
  using dsp      = butterworth<8>;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f = 500.;

  dsp::highpass ({test.scalar_coeffs}, f, test.samplerate);
  dsp::highpass_simd (
    test.vector_coeffs, vec_set<svf_test::vector_type> (f), test.samplerate);

  test();
}
//------------------------------------------------------------------------------
} // namespace artv
