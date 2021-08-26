#include <array>
#include <cstdio>
#include <cstring>
#include <gtest/gtest.h>

#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "common/simd_block_test.hpp"

namespace artv {

//------------------------------------------------------------------------------
TEST (svf, lowpass)
{
  using dsp      = andy::svf;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f = 500.;
  svf_test::scalar_type q = 2.;

  dsp::lowpass ({test.scalar_coeffs}, f, q, test.samplerate);
  dsp::lowpass_simd (
    test.vector_coeffs,
    vec_set<svf_test::vector_type> (f),
    vec_set<svf_test::vector_type> (q),
    test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (svf, highpass)
{
  using dsp      = andy::svf;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f = 500.;
  svf_test::scalar_type q = 2.;

  dsp::highpass ({test.scalar_coeffs}, f, q, test.samplerate);
  dsp::highpass_simd (
    test.vector_coeffs,
    vec_set<svf_test::vector_type> (f),
    vec_set<svf_test::vector_type> (q),
    test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (svf, allpass)
{
  using dsp      = andy::svf;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f = 500.;
  svf_test::scalar_type q = 2.;

  dsp::allpass ({test.scalar_coeffs}, f, q, test.samplerate);
  dsp::allpass_simd (
    test.vector_coeffs,
    vec_set<svf_test::vector_type> (f),
    vec_set<svf_test::vector_type> (q),
    test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (svf, bell)
{
  using dsp      = andy::svf;
  using svf_test = simd_block_test<dsp>;

  svf_test              test;
  svf_test::scalar_type f  = 500.;
  svf_test::scalar_type q  = 2.;
  svf_test::scalar_type db = 12.;

  andy::svf::bell ({test.scalar_coeffs}, f, q, db, test.samplerate);
  andy::svf::bell_simd (
    test.vector_coeffs,
    vec_set<svf_test::vector_type> (f),
    vec_set<svf_test::vector_type> (q),
    vec_set<svf_test::vector_type> (db),
    test.samplerate);

  test();
}
//------------------------------------------------------------------------------
} // namespace artv
