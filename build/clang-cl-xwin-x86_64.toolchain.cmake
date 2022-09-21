
# sudo apt install lld-15 clang-tools-15
# cmake -DCMAKE_TOOLCHAIN_FILE=./clang-cl-xwin-x86_64.toolchain.cmake -DCMAKE_BUILD_TYPE=Release -B win-release-build -S .

# Heavy lifting done here, just dumbing down:
# https://github.com/Nemirtingas/clang-msvc-sdk/blob/master/clang-cl-msvc.cmake

set (CMAKE_SYSTEM_NAME Windows)
set (CMAKE_SYSTEM_PROCESSOR x86_64)

set (CLANG_VERSION 15)

# cross tools
set (CMAKE_C_COMPILER clang-cl-${CLANG_VERSION})
set (CMAKE_CXX_COMPILER clang-cl-${CLANG_VERSION})
set (CMAKE_RC_COMPILER llvm-rc-${CLANG_VERSION})
set (CMAKE_AR llvm-lib-${CLANG_VERSION})
set (CMAKE_RANLIB llvm-ranlib-${CLANG_VERSION})
SET (CMAKE_MT "llvm-mt-${CLANG_VERSION}")
set (CMAKE_LINKER "lld-link-${CLANG_VERSION}")

# target environment on the build host system
set (CMAKE_FIND_ROOT_PATH /usr/bin)

# modify default behavior of FIND_XXX() commands
set (CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set (CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set (CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

set (XWIN_ROOT "${CMAKE_SOURCE_DIR}/xwin-0.2.5/files")
 # Compile flags ----
set (
COMPILE_FLAGS
    -D_CRT_SECURE_NO_WARNINGS
    -fuse-ld=lld
    --target=x86_64-windows-msvc
    # https://dev.to/yumetodo/list-of-mscver-and-mscfullver-8nd
    -fms-compatibility-version=19.29 # 2019 Update 11
    -Xclang -fexceptions
    -Xclang -fcxx-exceptions
    /winsdkdir "${XWIN_ROOT}/sdk"
    /vctoolsdir "${XWIN_ROOT}/crt"
    )
string (REPLACE ";" " " COMPILE_FLAGS "${COMPILE_FLAGS}")

# CMake rigmarole to preserve user-passed flags.
set(_CMAKE_C_FLAGS_INITIAL "${CMAKE_C_FLAGS}" CACHE STRING "")
set(
    CMAKE_C_FLAGS
    "${_CMAKE_C_FLAGS_INITIAL} ${COMPILE_FLAGS}"
    CACHE STRING "" FORCE
    )
# CMake rigmarole to preserve user-passed flags.
set(_CMAKE_CXX_FLAGS_INITIAL "${CMAKE_CXX_FLAGS}" CACHE STRING "")
set(
    CMAKE_CXX_FLAGS
    "${_CMAKE_CXX_FLAGS_INITIAL} ${COMPILE_FLAGS}"
    CACHE STRING "" FORCE
    )

# Link flags ----
set(LINK_FLAGS
    -libpath:${XWIN_ROOT}/crt/lib/x86_64
    -libpath:${XWIN_ROOT}/sdk/lib/um/x86_64
    -libpath:${XWIN_ROOT}/sdk/lib/ucrt/x86_64
    )
string (REPLACE ";" " " LINK_FLAGS "${LINK_FLAGS}")

# CMake rigmarole to preserve user-passed flags.
set(_CMAKE_EXE_LINKER_FLAGS_INITIAL "${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS "${_CMAKE_EXE_LINKER_FLAGS_INITIAL} ${LINK_FLAGS}" CACHE STRING "" FORCE)

# CMake rigmarole to preserve user-passed flags.
set(_CMAKE_MODULE_LINKER_FLAGS_INITIAL "${CMAKE_MODULE_LINKER_FLAGS}" CACHE STRING "")
set(CMAKE_MODULE_LINKER_FLAGS "${_CMAKE_MODULE_LINKER_FLAGS_INITIAL} ${LINK_FLAGS}" CACHE STRING "" FORCE)

# CMake rigmarole to preserve user-passed flags.
set(_CMAKE_SHARED_LINKER_FLAGS_INITIAL "${CMAKE_SHARED_LINKER_FLAGS}" CACHE STRING "")
set(CMAKE_SHARED_LINKER_FLAGS "${_CMAKE_SHARED_LINKER_FLAGS_INITIAL} ${LINK_FLAGS}" CACHE STRING "" FORCE)

# RC flags ----
set(_CMAKE_RC_FLAGS_INITIAL
    -I${XWIN_ROOT}/crt/include
    -I${XWIN_ROOT}/sdk/include/ucrt
    -I${XWIN_ROOT}/sdk/include/um
    -I${XWIN_ROOT}/sdk/include/shared)
string(REPLACE ";" " " _CMAKE_RC_FLAGS_INITIAL "${_CMAKE_RC_FLAGS_INITIAL}")
set(CMAKE_RC_FLAGS "${_CMAKE_RC_FLAGS_INITIAL}" CACHE STRING "" FORCE)
