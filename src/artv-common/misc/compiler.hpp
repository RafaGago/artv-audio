#pragma once

#include <bl/base/compiler.h>

#include "artv-common/misc/short_ints.hpp"

namespace artv {

//------------------------------------------------------------------------------
#if BL_COMPILER == BL_COMPILER_CLANG //-----------------------------------------
//------------------------------------------------------------------------------
template <uint bytes, class T>
auto assume_aligned_hint (T* ptr __attribute__ ((align_value (bytes))))
{
  return ptr;
}
//------------------------------------------------------------------------------
#define ARTV_LOOP_UNROLL_SIZE_HINT(n_loops) \
  ARTV_PRAGMA (clang loop unroll_count (n_loops))
//------------------------------------------------------------------------------
#define ARTV_PRAGMA(...) _Pragma (#__VA_ARGS__)
//------------------------------------------------------------------------------
#define artv_restrict __restrict__
//------------------------------------------------------------------------------
#elif BL_COMPILER == BL_COMPILER_GCC //-----------------------------------------
//------------------------------------------------------------------------------
template <uint bytes, class T>
auto assume_aligned_hint (T* ptr)
{
  return (T*) __builtin_aligned (ptr, bytes);
}
//------------------------------------------------------------------------------
#define ARTV_LOOP_UNROLL_SIZE_HINT(n_loops) ARTV_PRAGMA (GCC unroll n_loops)
//------------------------------------------------------------------------------
#define ARTV_PRAGMA(...) _Pragma (#__VA_ARGS__)
//------------------------------------------------------------------------------
#define artv_restrict __restrict__
//------------------------------------------------------------------------------
#else //------------------------------------------------------------------------
//------------------------------------------------------------------------------
T* assume_aligned_hint (T* ptr)
{
  // todo: warning?
  return ptr;
}
//------------------------------------------------------------------------------
#define ARTV_LOOP_UNROLL_SIZE_HINT(n_loops) ARTV_PRAGMA (GCC unroll n_loops)
//------------------------------------------------------------------------------
#define ARTV_PRAGMA(...)
//------------------------------------------------------------------------------
#define artv_restrict
//------------------------------------------------------------------------------
#endif

} // namespace artv