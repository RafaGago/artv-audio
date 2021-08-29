#pragma once

namespace artv {

//------------------------------------------------------------------------------
// https://www.musicdsp.org/en/latest/Synthesis/216-fast-whitenoise-generator.html
class white_noise_generator {
public:
  float operator() (float flevel) // will generate from -flevel to +flevel
  {
    flevel *= g_fScale;
    g_x1 ^= g_x2;
    float ret = g_x2 * flevel;
    g_x2 += g_x1;
    return ret;
  }
  // TODO: make a double version

private:
  static constexpr double g_fScale = 2.0f / 4294967296.f;

  int g_x1 = 0x67452301;
  int g_x2 = 0xefcdab89;
};
//------------------------------------------------------------------------------

} // namespace artv
