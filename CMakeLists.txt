# HINT for a recent version on Debian "pip3 install cmake"
cmake_minimum_required (VERSION 3.12.0)

include (FetchContent)
include (ExternalProject)
#-------------------------------------------------------------------------------
option (BUILD_VST3 "build VST3" ON)
option (BUILD_LV2 "build LV2" OFF) # Requires lv2-dev
option (LINK_TIME_OPTIMIZATION "Enable LTO" OFF) # To save build time turn off
# TODO: this one should be taken from the command line, but it wasn't
# trivial with clang-cl and the toolchain file.
option (ZEN3_BUILD "Build for the developer's DAW machine" OFF)
#-------------------------------------------------------------------------------
set (PLUGIN_FORMATS "")
if (BUILD_LV2)
    if (UNIX AND NOT APPLE)
        list (APPEND PLUGIN_FORMATS LV2)
    else()
        error ("LV2 Unsupported on this platform. Sorry.")
    endif()
endif()
if (BUILD_VST3)
    list (APPEND PLUGIN_FORMATS VST3)
endif()
#-------------------------------------------------------------------------------
function (configure_address_sanitizer target)
    # Linux only as of now.
    set (
        ASAN_FLAGS
        -fsanitize-address-use-after-scope
        -fsanitize-address-use-after-return=always
        -fsanitize=address
    )
    target_compile_options(
        ${target} PRIVATE -fno-omit-frame-pointer ${ASAN_FLAGS}
        )
    target_link_options (${target} PUBLIC ${ASAN_FLAGS})
endfunction()

function (configure_memory_sanitizer target)
    # As of now there are so many JUCE relatesd errors. It is not practical
    # to use.
    set (
        MEMSAN_FLAGS
        -fsanitize=memory
        -fsanitize-recover=all
        #-fsanitize-address-use-after-return=always
        -fsanitize-ignorelist=${CMAKE_SOURCE_DIR}/build/memory-sanitizer-ignorelist.txt
    )
    # Lots of positives on JUCE. Using ignorelist.
    target_compile_options(
        ${target} PUBLIC -fno-omit-frame-pointer ${MEMSAN_FLAGS}
        )
    target_link_options(
        ${target} PUBLIC ${MEMSAN_FLAGS}
        )
endfunction()
#-------------------------------------------------------------------------------
# This version number is mostly unused now...
set (VERSION_MAJOR 1)
set (VERSION_MINOR 0)
set (VERSION_REV 0)
set (VERSION_FULL "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_REV}")
set (CMAKE_CXX_STANDARD 17)

project (artv-audio VERSION ${VERSION_FULL})
set (INSTALL_DIR ${CMAKE_BINARY_DIR}/install)
set (VST3_COPY_DIR vst3-${CMAKE_BUILD_TYPE})

#-------------------------------------------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS On)
#-------------------------------------------------------------------------------
if (NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
    NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Notice that only "clang(cl)" and "gcc" are supported. Using GNU-style
    # flags on both.
    #
    # Parts of this codebase doesn't even compile on MSVC because of the
    # non-standards-compliant Visual Studio C preprocessor. It seems that MS
    # is doing an effort to bring its PP to compliance, so this codebase might
    # work on MSVC on the future.
    message (FATAL_ERROR "Unsuported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()
#-------------------------------------------------------------------------------
if ("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    # These compiler flags are set WITHOUT profiling, just a starting point
    # based on informed guesses.
    if (LINK_TIME_OPTIMIZATION)
        set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE) # LTO
        set (CMAKE_POLICY_DEFAULT_CMP0069 NEW # Always enable LTO (IPO, above)
        )
    endif()

    # Forced sse4.2 as of now. On the future maybe enabling AVX and other simd
    # sets could be done.
    # REMINDER: "-fopt-info-vec-missed" to see missed vectorizations.
    set (ARTV_DSP_COMPILE_OPTIONS_LIST -ffast-math -ftree-vectorize)

    set (ARTV_PLAIN_COMPILE_OPTIONS_LIST -Oz)
    # TODO: detecting target architecture is nontrivial from CMake, if compiling
    # for other arch than x64 the sse flags will need to be conditional
    set (ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST -O3 -DNDEBUG -msse4.2)

    if (ZEN3_BUILD)
        # TODO: this one should be taken from the command line, but it wasn't
        # trivial with clang-cl and the toolchain file.
        set (ZEN3_FLAGS -mtune=znver3 -march=znver3)
        set (ARTV_DSP_COMPILE_OPTIONS_LIST ${ARTV_DSP_COMPILE_OPTIONS_LIST} ${ZEN3_FLAGS})
        set (ARTV_PLAIN_COMPILE_OPTIONS_LIST ${ARTV_PLAIN_COMPILE_OPTIONS_LIST} ${ZEN3_FLAGS})
        set (ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST ${ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST} ${ZEN3_FLAGS})
    endif()
else()
    # Warnings on debug/dev builds
    set (ARTV_DSP_COMPILE_OPTIONS_LIST -ffast-math)
    set (ARTV_PLAIN_COMPILE_OPTIONS_LIST -ffast-math)
    set (ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST ${ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST})

    if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        # TODO: test CMAKE_MSVC_DEBUG_INFORMATION_FORMAT
        set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Od /MDd /Z7")
    else()
        set(
            CMAKE_EXE_LINKER_FLAGS
            ${CMAKE_EXE_LINKER_FLAGS}
            -Wl,--compress-debug-sections=zlib
            )
        set (
            CMAKE_MODULE_LINKER_FLAGS
            ${CMAKE_MODULE_LINKER_FLAGS}
            -Wl,--compress-debug-sections=zlib
            )
        set(
            CMAKE_SHARED_LINKER_FLAGS
            ${CMAKE_SHARED_LINKER_FLAGS}
            -Wl,--compress-debug-sections=zlib
            )
    endif()
endif()

set(HANA_CXX_FLAGS "-Wno-unused-value")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set (
        ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST
        ${ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST}
        -fdenormal-fp-math=positive-zero
        #-Rpass-analysis=loop-vectorize
        )
    if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        # convert arguments
        list (TRANSFORM ARTV_DSP_COMPILE_OPTIONS_LIST PREPEND "-clang:")
        list (TRANSFORM ARTV_PLAIN_COMPILE_OPTIONS_LIST PREPEND "-clang:")
        list (TRANSFORM ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST PREPEND "-clang:")
        list (TRANSFORM HANA_CXX_FLAGS PREPEND "-clang:")

        if ("${CMAKE_BUILD_TYPE}" MATCHES "Release")
            # TODO: test CMAKE_MSVC_DEBUG_INFORMATION_FORMAT
            # generate pdb
            set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi")
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /DEBUG /OPT:REF /OPT:ICF")
            set (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /DEBUG /OPT:REF /OPT:ICF")
            set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /DEBUG /OPT:REF /OPT:ICF")

            set (ARTV_DSP_COMPILE_OPTIONS_LIST ${ARTV_DSP_COMPILE_OPTIONS_LIST} /Zi)
            set (ARTV_PLAIN_COMPILE_OPTIONS_LIST ${ARTV_PLAIN_COMPILE_OPTIONS_LIST} /Zi)
            set (ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST ${ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST} /Zi)
        endif()
    else()
        # Always add debug symbols
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")

        set (ARTV_DSP_COMPILE_OPTIONS_LIST ${ARTV_DSP_COMPILE_OPTIONS_LIST} -g)
        set (ARTV_PLAIN_COMPILE_OPTIONS_LIST ${ARTV_PLAIN_COMPILE_OPTIONS_LIST} -g)
        set (ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST ${ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST} -g)
    endif()
endif()

set (CMAKE_CXX_FLAGS_RELEASE "${ARTV_CMAKE_CXX_FLAGS_RELEASE}")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftemplate-backtrace-limit=0")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")
set (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fuse-ld=lld")
set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=lld")

# message(FATAL_ERROR ${ARTV_DSP_COMPILE_OPTIONS_LIST})

# These targets are flag overrides that have to be added last, just to avoid
# dealing with CMake configuration messes which I refuse to learn and debug:
# https://stackoverflow.com/questions/24460486/cmake-build-type-is-not-being-used-in-cmakelists-txt/24470998#24470998
add_library (plain_target_compile_options INTERFACE)
target_compile_options(
    plain_target_compile_options INTERFACE ${ARTV_PLAIN_COMPILE_OPTIONS_LIST}
   )
add_library (dsp_target_compile_options INTERFACE)
target_compile_options(
    dsp_target_compile_options INTERFACE ${ARTV_DSP_COMPILE_OPTIONS_LIST}
    )

# CMake string/list awfulness . Depressing...
string(
    REPLACE
    ";" " "
    ARTV_DSP_COMPILE_OPTIONS
    "${ARTV_DSP_COMPILE_OPTIONS_LIST}"
    )
string(
    REPLACE
    ";" " "
    ARTV_PLAIN_COMPILE_OPTIONS
    "${ARTV_PLAIN_COMPILE_OPTIONS_LIST}"
    )
string(
    REPLACE
    ";" " "
    CMAKE_CXX_FLAGS_RELEASE
    "${ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST}"
    )

#-- Juce -----------------------------------------------------------------------
set (JUCE_PATCH ${CMAKE_CURRENT_SOURCE_DIR}/patches/juce/0001-Slider-Add-fine-scrolling.patch)

add_compile_definitions(
    JUCER_ENABLE_GPL_MODE=1
    JUCE_DISPLAY_SPLASH_SCREEN=0
    JUCE_USE_SSE_INTRINSICS=1
    JUCE_DISABLE_CAUTIOUS_PARAMETER_ID_CHECKING=1 # To disable AAX-related asserts
    )
# Juce 7.0.5
FetchContent_Declare(
    juce
    GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
    GIT_TAG        69795dc8e589a9eb5df251b6dd994859bf7b3fab
    PATCH_COMMAND  git reset --hard 965d0ca4be178c4a0000b116d460e15c30311992 && git apply --whitespace=fix ${JUCE_PATCH}
    GIT_SHALLOW    ON
    )
FetchContent_MakeAvailable (juce)
set (JUCE_REPO ${juce_SOURCE_DIR})

#-- Fonts ----------------------------------------------------------------------
add_subdirectory(resources/fonts)

#-- Cockos WDL -----------------------------------------------------------------
FetchContent_Declare(
    wdl_fetch
    GIT_REPOSITORY https://github.com/justinfrankel/WDL.git
    GIT_TAG        3cd6f4d11b5986e7fb6975a66c96f9cbef5f5f37
    )
FetchContent_MakeAvailable (wdl_fetch)
set (WDL_REPO ${wdl_fetch_SOURCE_DIR})

add_library (wdl INTERFACE)
target_include_directories (wdl INTERFACE ${WDL_REPO})

#-- Own Utilities library (C) --------------------------------------------------

# Just using the headers, avoiding meson for now.

set (BL_REPO ${CMAKE_BINARY_DIR}/bl-git)
ExternalProject_Add(
    bl_fetch
    GIT_REPOSITORY https://github.com/RafaGago/base_library.git
    GIT_TAG        3d08c1ccbecc4ab26495999fc6692978b490e40f
    #GIT_TAG        7.6.0
    SOURCE_DIR     "${BL_REPO}"
    UPDATE_COMMAND    ""
    BUILD_ALWAYS   False
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )
add_library (bl_headers INTERFACE)
target_include_directories (bl_headers INTERFACE ${BL_REPO}/include)
add_dependencies (bl_headers bl_fetch)
#-- Boost hanna-----------------------------------------------------------------
set (HANA_REPO ${CMAKE_BINARY_DIR}/boost-hana-git)

ExternalProject_Add(
    hana_fetch
    GIT_REPOSITORY https://github.com/boostorg/hana.git
    GIT_TAG        v1.7.0
    SOURCE_DIR     "${HANA_REPO}"
    UPDATE_COMMAND    ""
    BUILD_ALWAYS   False
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )
add_library (hana INTERFACE)
target_include_directories (hana INTERFACE ${HANA_REPO}/include)
add_dependencies (hana hana_fetch)
# TODO check if future versions remove the need for this supression
target_compile_options(hana INTERFACE ${HANA_CXX_FLAGS})
#-- Boost mp11------------------------------------------------------------------
set (MP11_REPO ${CMAKE_BINARY_DIR}/boost-mp11-git)

ExternalProject_Add(
    mp11_fetch
    GIT_REPOSITORY https://github.com/boostorg/mp11.git
    GIT_TAG        boost-1.77.0
    SOURCE_DIR     "${MP11_REPO}"
    UPDATE_COMMAND    ""
    BUILD_ALWAYS   False
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )

add_library (mp11 INTERFACE)
target_include_directories (mp11 INTERFACE ${MP11_REPO}/include)
add_dependencies (mp11 mp11_fetch)
#-- XSIMD ----------------------------------------------------------------------
set (XSIMD_REPO ${CMAKE_BINARY_DIR}/xsimd-git)
ExternalProject_Add(
    xsimd_fetch
    GIT_REPOSITORY https://github.com/RafaGago/xsimd.git
    GIT_TAG        1075472af267695730adba03311d5b62b266cdf9
    SOURCE_DIR     "${XSIMD_REPO}"
    UPDATE_COMMAND    ""
    BUILD_ALWAYS   False
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )
add_library (xsimd INTERFACE)
target_include_directories (xsimd INTERFACE ${XSIMD_REPO}/include)
add_dependencies (xsimd xsimd_fetch)

#-- Gtest ----------------------------------------------------------------------
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.10.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable (googletest)
#-- Gbenchmark -----------------------------------------------------------------
# FetchContent_Declare(
#     googlebenchmark
#     GIT_REPOSITORY https://github.com/google/benchmark.git
#     GIT_TAG        v1.5.0
# )
# set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
# set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable (googlebenchmark)
#--muFFT (float)----------------------------------------------------------------
set (MUFFT_INCLUDE ${CMAKE_BINARY_DIR}/mufft-src-hack)
set (MUFFT_REPO ${MUFFT_INCLUDE}/mufft)

FetchContent_Declare(
    mufft_fetch
    GIT_REPOSITORY https://github.com/Themaister/muFFT.git
    GIT_TAG        47bb08652eab399c2c7d460abe5184857110f130
    SOURCE_DIR     "${MUFFT_REPO}"
    )
FetchContent_MakeAvailable (mufft_fetch)

set_target_properties (muFFT PROPERTIES POSITION_INDEPENDENT_CODE On)
set_target_properties (muFFT-sse PROPERTIES POSITION_INDEPENDENT_CODE On)
set_target_properties (muFFT-sse3 PROPERTIES POSITION_INDEPENDENT_CODE On)
set_target_properties (muFFT-avx PROPERTIES POSITION_INDEPENDENT_CODE On)

add_library (mufft-all INTERFACE)
target_include_directories (mufft-all INTERFACE ${MUFFT_INCLUDE})
target_link_libraries (mufft-all INTERFACE muFFT muFFT-sse muFFT-sse3 muFFT-avx)
#--PFFFT port (double) ----------------------------------------------------------

set (PFFFT_INCLUDE ${CMAKE_BINARY_DIR}/pffft-git)
set (PFFFT_REPO ${PFFFT_INCLUDE}/pffft) # to be able to "#include" cleanly.
ExternalProject_Add(
    pffft_fetch
    GIT_REPOSITORY    https://github.com/marton78/pffft.git
    GIT_TAG           9603871ac44b6f782e9666a968d25977f7124e39
    GIT_SUBMODULES    ""
    SOURCE_DIR        "${PFFFT_REPO}"
    UPDATE_COMMAND    ""
    BUILD_ALWAYS      False
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )

# Ignoring the provided CMakeLists with hardcoded flags, just going after the C
# files.

add_library (pffft STATIC)
target_compile_features (pffft PUBLIC c_std_99)
target_include_directories (pffft PUBLIC ${PFFFT_INCLUDE})
target_include_directories (pffft PRIVATE ${PFFFT_REPO})
target_sources(
    pffft
    PRIVATE
    ${PFFFT_REPO}/pffft_double.c
    ${PFFFT_REPO}/pffft_common.c
    )
set_source_files_properties(
    ${PFFFT_REPO}/pffft_double.c
    ${PFFFT_REPO}/pffft_common.c
    PROPERTIES GENERATED 1
    )
target_compile_definitions (pffft PRIVATE _USE_MATH_DEFINES)

set_property (TARGET pffft PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries (pffft PRIVATE dsp_target_compile_options)
add_dependencies (pffft pffft_fetch)

#--GCEM (constexpr math) -------------------------------------------------------
set (GCEM_REPO ${CMAKE_BINARY_DIR}/gcem)

ExternalProject_Add(
    gcem_fetch
    GIT_REPOSITORY https://github.com/kthohr/gcem.git
    GIT_TAG        a20b0fc0206ff7d99a96fe4afdfe8205b01c8220
    SOURCE_DIR     "${GCEM_REPO}"
    UPDATE_COMMAND    ""
    BUILD_ALWAYS   False
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )
add_library (gcem INTERFACE)
target_include_directories (gcem INTERFACE ${GCEM_REPO}/include)
add_dependencies (gcem gcem_fetch)

#-- DISHTRO Ports  -------------------------------------------------------------
set (DISTHRO_PORTS_INCLUDE ${CMAKE_BINARY_DIR}/disthro-ports-root)
# Adding one intermediate directory.
set (DISTHRO_PORTS_REPO ${DISTHRO_PORTS_INCLUDE}/disthro-ports)
ExternalProject_Add(
    dishthro_ports_fetch
    GIT_REPOSITORY https://github.com/RafaGago/DISTRHO-Ports.git
    GIT_TAG        6248eda3f490e30a0f4c13e6b0a04dfa120022cf
    SOURCE_DIR     "${DISTHRO_PORTS_REPO}"
    UPDATE_COMMAND    ""
    BUILD_ALWAYS   False
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )
add_library (disthro_ports INTERFACE)
target_include_directories (disthro_ports INTERFACE ${DISTHRO_PORTS_INCLUDE})
add_dependencies (disthro_ports dishthro_ports_fetch)

add_library (tal_reverb2 INTERFACE)
target_link_libraries (tal_reverb2 INTERFACE disthro_ports)

add_library (luftikus INTERFACE)
target_sources(
    luftikus
    INTERFACE
    ${DISTHRO_PORTS_REPO}/ports-legacy/luftikus/source/dsp/eqdsp.cpp
    )
set_source_files_properties(
    ${DISTHRO_PORTS_REPO}/ports-legacy/luftikus/source/dsp/eqdsp.cpp
    PROPERTIES
    GENERATED 1
    )
target_link_libraries (luftikus INTERFACE disthro_ports)

# -- DragonFly reverb #----------------------------------------------------------
FetchContent_Declare(
    dragonfly_fetch
    GIT_REPOSITORY https://github.com/RafaGago/dragonfly-reverb.git
    GIT_TAG        dsp-namespaces
    )
FetchContent_MakeAvailable (dragonfly_fetch)
set (DRAGONFLY_REPO ${dragonfly_fetch_SOURCE_DIR})

set (DRAGONFLY_FREEVERB_PATH ${DRAGONFLY_REPO}/common/freeverb)
add_library (dragonfly_freeverb STATIC)
target_include_directories (
    dragonfly_freeverb PRIVATE ${DRAGONFLY_FREEVERB_PATH}/..
    )
target_sources(
    dragonfly_freeverb
    PRIVATE
    ${DRAGONFLY_FREEVERB_PATH}/allpass.cpp
    ${DRAGONFLY_FREEVERB_PATH}/biquad.cpp
    ${DRAGONFLY_FREEVERB_PATH}/comb.cpp
    ${DRAGONFLY_FREEVERB_PATH}/delay.cpp
    ${DRAGONFLY_FREEVERB_PATH}/delayline.cpp
    ${DRAGONFLY_FREEVERB_PATH}/earlyref.cpp
    ${DRAGONFLY_FREEVERB_PATH}/efilter.cpp
    ${DRAGONFLY_FREEVERB_PATH}/nrev.cpp
    ${DRAGONFLY_FREEVERB_PATH}/nrevb.cpp
    ${DRAGONFLY_FREEVERB_PATH}/progenitor.cpp
    ${DRAGONFLY_FREEVERB_PATH}/progenitor2.cpp
    ${DRAGONFLY_FREEVERB_PATH}/revbase.cpp
    ${DRAGONFLY_FREEVERB_PATH}/slot.cpp
    ${DRAGONFLY_FREEVERB_PATH}/strev.cpp
    ${DRAGONFLY_FREEVERB_PATH}/utils.cpp
    ${DRAGONFLY_FREEVERB_PATH}/zrev.cpp
    ${DRAGONFLY_FREEVERB_PATH}/zrev2.cpp
    )
set_property (TARGET dragonfly_freeverb PROPERTY CXX_STANDARD 98)
set_property (TARGET dragonfly_freeverb PROPERTY POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(
    dragonfly_freeverb PRIVATE LIBFV3_FLOAT=1
    )
target_link_libraries (dragonfly_freeverb PRIVATE dsp_target_compile_options)

function (add_dragonfly_plugin_lib target plugin_name)
    add_library (${target} STATIC)
    set_property (TARGET ${target} PROPERTY CXX_STANDARD 14)
    set_property (TARGET ${target} PROPERTY POSITION_INDEPENDENT_CODE ON)

    set(DRAGONFLY_PLATFORM_DEFINITIONS "")
    if(MSVC)
        # MSVC Includes clang-cl. Used because of Dragonfly/Dishtro
        set(DRAGONFLY_PLATFORM_DEFINITIONS "ssize_t=uint64_t")
    endif()

    target_compile_definitions(
        ${target}
        PRIVATE
        LIBFV3_FLOAT=1
        HAVE_CPP11_SUPPORT=1
        ${DRAGONFLY_PLATFORM_DEFINITIONS}
        )
    target_sources(
        ${target}
        PRIVATE
        ${DRAGONFLY_REPO}/plugins/${plugin_name}/DSP.cpp
        src/artv-common/dsp/third_party/dragonfly/${target}_cfirewall.cpp
        )
    target_include_directories(
        ${target}
        PRIVATE
        ${DRAGONFLY_REPO}/common
        ${DRAGONFLY_REPO}/dpf/distrho
        ${DRAGONFLY_REPO}/plugins
        ${DRAGONFLY_REPO}/plugins/${plugin_name}
        src
        )
    target_link_libraries (${target} PRIVATE dsp_target_compile_options)
endfunction()

add_dragonfly_plugin_lib (dragonfly_er    dragonfly-early-reflections)
add_dragonfly_plugin_lib (dragonfly_plate dragonfly-plate-reverb)
add_dragonfly_plugin_lib (dragonfly_hall  dragonfly-hall-reverb)
add_dragonfly_plugin_lib (dragonfly_room  dragonfly-room-reverb)

add_library (dragonfly INTERFACE)
target_link_libraries(
    dragonfly
    INTERFACE
    dragonfly_er
    dragonfly_plate
    dragonfly_room
    dragonfly_hall
    dragonfly_freeverb
    )
#-- FF-meters-------------------------------------------------------------------
FetchContent_Declare(
    ff_meters_fetch
    GIT_REPOSITORY https://github.com/RafaGago/ff_meters.git
    GIT_TAG        6ccdab82ca521115dddabb351a9ced7e4daab35a
    )
FetchContent_MakeAvailable (ff_meters_fetch)
set (FF_METERS_REPO ${ff_meters_fetch_SOURCE_DIR})

add_library (ff_meters INTERFACE)
target_include_directories (ff_meters INTERFACE ${FF_METERS_REPO})
target_compile_definitions(
    ff_meters INTERFACE FF_AUDIO_ALLOW_ALLOCATIONS_IN_MEASURE_BLOCK=0)
# Don't want to deal on how to compile a static lib based on JUCE, so just
# adding the sources to the relevant project via INTERFACE.
target_sources(
    ff_meters
    INTERFACE
    ${FF_METERS_REPO}/ff_meters.cpp
   )
set_source_files_properties(
    ${FF_METERS_REPO}/ff_meters.cpp
    PROPERTIES
    COMPILE_FLAGS
    "${ARTV_PLAIN_COMPILE_OPTIONS}"
    )
#-- Rubberband -----------------------------------------------------------------
FetchContent_Declare(
    rubberband_fetch
    GIT_REPOSITORY https://github.com/breakfastquay/rubberband
    GIT_TAG        41ba70fbb9dc5042da37395a9d7a28fd967036d0
    )
FetchContent_MakeAvailable (rubberband_fetch)
set (RUBBERBAND_REPO ${rubberband_fetch_SOURCE_DIR})

add_library (rubberband STATIC)
target_include_directories (rubberband INTERFACE ${RUBBERBAND_REPO})
target_sources(
    rubberband
    PRIVATE
    ${RUBBERBAND_REPO}/single/RubberBandSingle.cpp
   )
target_compile_definitions (rubberband PRIVATE NO_THREADING=1)
set_property (TARGET rubberband PROPERTY POSITION_INDEPENDENT_CODE ON)
# Avoid some problems with ambiguous references to "byte" by setting the
# standard to C++14.
set_property (TARGET rubberband PROPERTY CXX_STANDARD 14)
target_compile_definitions (rubberband PRIVATE NOMINMAX)
target_link_libraries (rubberband PRIVATE dsp_target_compile_options)

#-- Soundtouch -----------------------------------------------------------------
set (SOUNDTOUCH_PATCH ${CMAKE_CURRENT_SOURCE_DIR}/patches/soundtouch/0001-Remove-problematic-flags.patch)
FetchContent_Declare(
    soundtouch_fetch
    GIT_REPOSITORY https://codeberg.org/soundtouch/soundtouch.git
    GIT_TAG        9fedba866e36edfca11e51c852d92b903f0735aa
    PATCH_COMMAND  git reset --hard 9fedba866e36edfca11e51c852d92b903f0735aa && git apply --whitespace=fix ${SOUNDTOUCH_PATCH}
    )
FetchContent_MakeAvailable (soundtouch_fetch)
set_property (TARGET SoundTouch PROPERTY POSITION_INDEPENDENT_CODE ON)

# This is not done by linking against the "dsp_target_compile_options" interface
# library because that required bumping the CMAKE_REQUIRED_VERSION and to set
# the CMP0069 policy. This was not a problem locally, but it failed on GitHub's
# servers which are slow and painful to Debug (TODO if feeling veeeery bored).
target_compile_options(
    SoundTouch PRIVATE ${ARTV_DSP_COMPILE_OPTIONS_LIST}
    )
#-- Common ---------------------------------------------------------------------
add_library(common INTERFACE)
target_include_directories (common INTERFACE src)
target_link_libraries(
    common
    INTERFACE
    hana
    mp11
    xsimd
    gcem
    bl_headers
    mufft-all
    pffft)

if(MSVC)
    # MSVC Includes clang-cl.
    set(COMMON_PLATFORM_DEFINITIONS "_USE_MATH_DEFINES=1") # _M_PI, etc.
endif()
target_compile_definitions (common INTERFACE ${COMMON_PLATFORM_DEFINITIONS})

#-- Test- ----------------------------------------------------------------------
add_library (test INTERFACE)
target_include_directories (test INTERFACE test)
target_link_libraries(
    test
    INTERFACE
    gtest
    gtest_main
    common
    dsp_target_compile_options
    )
#-- Headless -------------------------------------------------------------------
add_library (headless INTERFACE)
target_include_directories (headless INTERFACE common)
target_include_directories (headless INTERFACE headless)
#-- Benchmark- -----------------------------------------------------------------
add_library (own_benchmark INTERFACE)
target_include_directories (own_benchmark INTERFACE benchmark)
target_link_libraries(
    own_benchmark
    INTERFACE
    benchmark
    benchmark_main
    common
    )
#-- Common tests ---------------------------------------------------------------
add_executable(
    common-test
    test/common/delay_compensation_buffers.cpp
    test/common/float.cpp
    test/common/fft.cpp
    test/common/simd.cpp
    test/common/fixed_point.cpp
    test/common/turbopaco.cpp
    )
target_link_libraries (common-test PUBLIC test mufft-all pffft)
#-- Mix-maxtrix ----------------------------------------------------------------
set (MIXMAXTRIX_VERSION_MAJOR 1)
set (MIXMAXTRIX_VERSION_MINOR 10)
set (MIXMAXTRIX_VERSION_REV 2)
set(
    MIXMAXTRIX_VERSION
    "${MIXMAXTRIX_VERSION_MAJOR}.${MIXMAXTRIX_VERSION_MINOR}.${MIXMAXTRIX_VERSION_REV}"
    )
function (add_mix_maxtrix_plugin_target target channels)
    juce_add_plugin(
        ${target}
        PLUGIN_MANUFACTURER_CODE ArtV
        PLUGIN_CODE Mmt${channels}
        PRODUCT_NAME ${target}
        LV2_SHARED_LIBRARY_NAME ${target}
        COMPANY_NAME ArtV
        VERSION "${MIXMAXTRIX_VERSION}"
        VST3_CATEGORIES Tools
        VST3_COPY_DIR ${VST3_COPY_DIR}
        FORMATS ${PLUGIN_FORMATS}
        NEEDS_MIDI_INPUT 0
        VST_NUM_MIDI_INS 0
        LV2_URI https://github.com/RafaGago/artv-audio
        )
    target_sources(
        ${target}
        PRIVATE
        src/artv-common/dsp/third_party/jsfx_engine/jsfx_engine.cpp
        src/mix-maxtrix/processor.cpp
        src/mix-maxtrix/editor.cpp
        )
    set_source_files_properties(
        src/mix-maxtrix/editor.cpp
        PROPERTIES
        COMPILE_FLAGS
        "${ARTV_PLAIN_COMPILE_OPTIONS}"
        )
    set_source_files_properties(
        src/artv-common/dsp/third_party/jsfx_engine/jsfx_engine.cpp
        src/mix-maxtrix/processor.cpp
        PROPERTIES
        COMPILE_FLAGS
        "${ARTV_DSP_COMPILE_OPTIONS}"
        )
    target_compile_definitions(
        ${target}
        PUBLIC
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
        JUCER_ENABLE_GPL_MODE=1
        VERSION_TXT="${MIXMAXTRIX_VERSION}"
        VERSION_MAJOR=${MIXMAXTRIX_VERSION_MAJOR}
        VERSION_MINOR=${MIXMAXTRIX_VERSION_MINOR}
        VERSION_REV=${MIXMAXTRIX_VERSION_REV}
        )
    target_link_libraries(
        ${target}
        PUBLIC
        common
        juce::juce_audio_utils
        dragonfly
        tal_reverb2
        luftikus
        ff_meters
        rubberband
        SoundTouch
        wdl
        )
    install(
        TARGETS ${target}
        RUNTIME DESTINATION plugins
        LIBRARY DESTINATION plugins
        ARCHIVE DESTINATION plugins
        )
endfunction()
add_mix_maxtrix_plugin_target (MixMaxTrix 8)
#-- Mix-maxtrix tests-----------------------------------------------------------
add_executable(
    mix-maxtrix-test
    test/mix-maxtrix/order_and_buffering.cpp
    )
target_link_libraries (mix-maxtrix-test PUBLIC test)
if(NOT MSVC)
    # This is in fact a GCC detection, good enough for now
    configure_address_sanitizer (mix-maxtrix-test)
endif()

add_executable(
    mix-maxtrix-smoke
    test/mix-maxtrix/smoke.cpp
    )
target_link_libraries (mix-maxtrix-smoke PUBLIC test MixMaxTrix)
if(NOT MSVC)
    # This is in fact a GCC detection, good enough for now
    configure_address_sanitizer (mix-maxtrix-smoke)
endif()
#-- Mix-maxtrix headless -------------------------------------------------------
# Necessary?
#juce_add_console_app(
#    MixMaxTrix-headless
#    PRODUCT_NAME MixMaxTrix
#    COMPANY_NAME Francisco
#    )
#juce_generate_juce_header (MixMaxTrix-headless)
#target_sources(
#    MixMaxTrix-headless
#    PRIVATE
#    headless/mix-maxtrix/main.cpp
#    )
#target_link_libraries (MixMaxTrix-headless PRIVATE headless MixMaxTrix)
#-- TurboPaco ----------------------------------------------------------------
set (TURBOPACO_VERSION_MAJOR 1)
set (TURBOPACO_VERSION_MINOR 0)
set (TURBOPACO_VERSION_REV 2)
set(
TURBOPACO_VERSION
"${TURBOPACO_VERSION_MAJOR}.${TURBOPACO_VERSION_MINOR}.${TURBOPACO_VERSION_REV}"
)
juce_add_plugin(
    turbopaco
    PLUGIN_MANUFACTURER_CODE ArtV
    PLUGIN_CODE Lfvb
    PRODUCT_NAME TurboPaco
    LV2_SHARED_LIBRARY_NAME turbopaco
    COMPANY_NAME ArtV
    VERSION "${TURBOPACO_VERSION}"
    VST3_CATEGORIES Reverb
    VST3_COPY_DIR ${VST3_COPY_DIR}
    FORMATS ${PLUGIN_FORMATS}
    NEEDS_MIDI_INPUT 0
    VST_NUM_MIDI_INS 0
    LV2_URI https://github.com/RafaGago/artv-audio
    )
target_sources(
    turbopaco
    PRIVATE
    src/turbopaco/processor.cpp
    src/turbopaco/editor.cpp
    )
set_source_files_properties(
    src/turbopaco/editor.cpp
    PROPERTIES
    COMPILE_FLAGS
    "${ARTV_PLAIN_COMPILE_OPTIONS}"
    )
set_source_files_properties(
    src/turbopaco/processor.cpp
    PROPERTIES
    COMPILE_FLAGS
    "${ARTV_DSP_COMPILE_OPTIONS} ${ARTV_TURBOPACO_COMPILE_OPTIONS}"
    )
target_compile_definitions(
    turbopaco
    PUBLIC
    JUCE_WEB_BROWSER=0
    JUCE_USE_CURL=0
    JUCE_VST3_CAN_REPLACE_VST2=0
    JUCER_ENABLE_GPL_MODE=1
    VERSION_TXT="${TURBOPACO_VERSION}"
    VERSION_MAJOR=${TURBOPACO_VERSION_MAJOR}
    VERSION_MINOR=${TURBOPACO_VERSION_MINOR}
    VERSION_REV=${TURBOPACO_VERSION_REV}
    )
target_link_libraries(
    turbopaco
    PUBLIC
    common
    juce::juce_audio_utils
    fonts
    )
install(
    TARGETS turbopaco
    RUNTIME DESTINATION plugins
    LIBRARY DESTINATION plugins
    ARCHIVE DESTINATION plugins
    )
