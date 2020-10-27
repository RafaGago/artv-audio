# HINT for a recent version on Debian "pip3 install cmake"
cmake_minimum_required(VERSION 3.12.0)

include (FetchContent)
include (ExternalProject)
#-------------------------------------------------------------------------------
option (BUILD_VST3 "build VST3" ON)
option (BUILD_LV2 "build LV2" OFF) # Requires lv2-dev
option (LINK_TIME_OPTIMIZATION "Enable LTO" OFF) # To save build time turn off

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
function (configure_memory_sanitizer target)
    # Linux only as of now.
    target_compile_options(
        ${target}
        PRIVATE
        -fsanitize=address -fno-omit-frame-pointer
        )
    target_link_options (${target} PUBLIC -fsanitize=address)
    target_compile_options(
        ${target}
        PUBLIC
        -fsanitize=address -fno-omit-frame-pointer
        )
    target_link_options(${target} PUBLIC -fsanitize=address)
endfunction()
#-------------------------------------------------------------------------------
set (VERSION_MAJOR 0)
set (VERSION_MINOR 0)
set (VERSION_REV 1)
set (VERSION_FULL "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_REV}")
set (CMAKE_CXX_STANDARD 17)

project (artv-audio VERSION ${VERSION_FULL})
set (INSTALL_DIR ${CMAKE_BINARY_DIR}/install)
set (VST3_COPY_DIR vst3-${CMAKE_BUILD_TYPE})
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

    # Forced sse2 as of now. On the future maybe enabling AVX and other simd
    # sets could be done.
    # REMINDER: "-fopt-info-vec-missed" to see missed vectorizations.
    set (ARTV_DSP_COMPILE_OPTIONS_LIST -ffast-math -ftree-vectorize)
    set (ARTV_PLAIN_COMPILE_OPTIONS_LIST -Os)
    set (ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST -O2 -DNDEBUG -msse2)
else()
    # Warnings on debug/dev builds
    set (ARTV_DSP_COMPILE_OPTIONS_LIST -ffast-math)
    set (ARTV_PLAIN_COMPILE_OPTIONS_LIST -ffast-math)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set (
        ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST
        ${ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST} -fdenormal-fp-math=positive-zero
        )
    if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        list (TRANSFORM ARTV_DSP_COMPILE_OPTIONS   PREPEND "/clang:")
        list (TRANSFORM ARTV_PLAIN_COMPILE_OPTIONS PREPEND "/clang:")
        list (TRANSFORM ARTV_CMAKE_CXX_FLAGS_RELEASE_LIST PREPEND "/clang:")
    endif()
endif()

set (CMAKE_CXX_FLAGS_RELEASE ${ARTV_CMAKE_CXX_FLAGS_RELEASE})

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
    REPLACE ";" " " ARTV_DSP_COMPILE_OPTIONS "${ARTV_DSP_COMPILE_OPTIONS_LIST}"
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
add_compile_definitions(
    JUCER_ENABLE_GPL_MODE=1
    JUCE_DISPLAY_SPLASH_SCREEN=0
    JUCE_USE_SSE_INTRINSICS=1
    )
FetchContent_Declare(
    juce
    GIT_REPOSITORY https://github.com/lv2-porting-project/JUCE
    GIT_TAG        f00a420bc348d79c4688b33b2b905b8ca0f25a3a
    )
FetchContent_MakeAvailable (juce)
set (JUCE_REPO ${juce_SOURCE_DIR})

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

# As this is a C library, make sure it inherits the release CMake flags for
# vectorization etc if used as main fft.
add_library (wdl-fft)
target_link_libraries (wdl-fft PUBLIC wdl)
target_link_libraries (wdl-fft PRIVATE dsp_target_compile_options)
target_sources (wdl-fft PRIVATE ${WDL_REPO}/WDL/fft.c)
target_include_directories (wdl-fft PRIVATE ${WDL_REPO}/WDL)
target_compile_definitions(
    wdl-fft
    PUBLIC
    WDL_FFT_REALSIZE=4
    )
#-- Boost hanna-----------------------------------------------------------------
set (HANA_REPO ${CMAKE_BINARY_DIR}/boost-hana-git)
set (HANA_INSTALL ${CMAKE_BINARY_DIR}/boost-hana-install)
ExternalProject_Add(
    hana_fetch
    GIT_REPOSITORY https://github.com/boostorg/hana.git
    GIT_TAG        v1.7.0
    SOURCE_DIR     "${HANA_REPO}"
    UPDATE_COMMAND ""
    BUILD_ALWAYS   False
    CMAKE_ARGS     "-DCMAKE_INSTALL_PREFIX=${HANA_INSTALL}"
    INSTALL_DIR    "${HANA_INSTALL}"
    )
add_library (hana INTERFACE)
target_include_directories (hana INTERFACE ${HANA_INSTALL}/include)
add_dependencies (hana hana_fetch)
#-- Boost mp11------------------------------------------------------------------
set (MP11_REPO ${CMAKE_BINARY_DIR}/boost-mp11-git)
set (MP11_INSTALL ${CMAKE_BINARY_DIR}/boost-mp11-install)
ExternalProject_Add(
    mp11_fetch
    GIT_REPOSITORY https://github.com/boostorg/mp11.git
    GIT_TAG        boost-1.75.0
    SOURCE_DIR     "${MP11_REPO}"
    UPDATE_COMMAND ""
    BUILD_ALWAYS   False
    CMAKE_ARGS     "-DCMAKE_INSTALL_PREFIX=${MP11_INSTALL}"
    INSTALL_DIR    "${MP11_INSTALL}"
    )
add_library (mp11 INTERFACE)
target_include_directories (mp11 INTERFACE ${MP11_INSTALL}/include)
add_dependencies (mp11 mp11_fetch)
#-- XSIMD ----------------------------------------------------------------------
set (XSIMD_REPO ${CMAKE_BINARY_DIR}/xsimd-git)
ExternalProject_Add(
    xsimd_fetch
    GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
    GIT_TAG        e8234f75b54c8af9f26b9f90a2ecdf5a22326fa9
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
#-- Hiir ----- -----------------------------------------------------------------
set (HIIR_REPO ${CMAKE_BINARY_DIR}/hiir-git)
ExternalProject_Add(
    hiir_fetch
    GIT_REPOSITORY https://github.com/unevens/hiir.git
    GIT_TAG        73a7b1f1c9a50f56fe85d71313223c5d4e61e889
    SOURCE_DIR     "${HIIR_REPO}"
    UPDATE_COMMAND    ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    )
add_library (hiir INTERFACE)
target_include_directories (hiir INTERFACE ${HIIR_REPO})
add_dependencies (hiir hiir_fetch)
#--muFFT -----------------------------------------------------------------------
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

#--GCEM (constexpr math) -------------------------------------------------------
set (GCEM_REPO ${CMAKE_BINARY_DIR}/gcem)

ExternalProject_Add(
    gcem_fetch
    GIT_REPOSITORY https://github.com/kthohr/gcem.git
    GIT_TAG        a20b0fc0206ff7d99a96fe4afdfe8205b01c8220
    SOURCE_DIR     "${GCEM_REPO}"
    UPDATE_COMMAND    ""
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
    GIT_TAG        3231d01751e98beaeca0f0d5d357cec5d6f50503
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
        src/artv-common/dsp/dragonfly/${target}_cfirewall.cpp
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
#-- Common ---------------------------------------------------------------------
add_library(common INTERFACE)
target_include_directories (common INTERFACE src)
target_link_libraries(common INTERFACE hana mp11 xsimd gcem)

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
    test/common/fft.cpp
    test/common/delay_compensation_buffers.cpp
    )
target_link_libraries (common-test PUBLIC test mufft-all wdl-fft)
#-- Mix-maxtrix ----------------------------------------------------------------
set (MIXMAXTRIX_VERSION_MAJOR 0)
set (MIXMAXTRIX_VERSION_MINOR 0)
set (MIXMAXTRIX_VERSION_REV 1)
set (MIXMAXTRIX_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_REV}")

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
        LV2_URI https://github.com/RafaGago/artv-audio-stuff
        )
    target_sources(
        ${target}
        PRIVATE
        src/artv-common/dsp/jsfx_engine/jsfx_engine.cpp
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
        src/artv-common/dsp/jsfx_engine/jsfx_engine.cpp
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
        MIXMAXTRIX_CHANNELS=${channels}
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
        mufft-all
        dragonfly
        tal_reverb2
        luftikus
        ff_meters
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
    configure_memory_sanitizer (mix-maxtrix-test)
endif()

add_executable(
    mix-maxtrix-smoke
    test/mix-maxtrix/smoke.cpp
    )
target_link_libraries (mix-maxtrix-smoke PUBLIC test MixMaxTrix)
if(NOT MSVC)
    # This is in fact a GCC detection, good enough for now
    configure_memory_sanitizer (mix-maxtrix-smoke)
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
