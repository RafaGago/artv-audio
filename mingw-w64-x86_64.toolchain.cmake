
# sudo apt install mingw-w64
# cmake -DCMAKE_TOOLCHAIN_FILE=./mingw-w64-x86_64.toolchain.cmake -DCMAKE_BUILD_TYPE=Release -B win-release-build -S .

set (CMAKE_SYSTEM_NAME Windows)
set (CMAKE_SYSTEM_PROCESSOR x86_64)
set (TOOLCHAIN_PREFIX x86_64-w64-mingw32)

# cross compilers to use for C, C++ and Fortran
set (CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc-posix)
set (CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++-posix)
set (CMAKE_RC_COMPILER ${TOOLCHAIN_PREFIX}-windres)
set (CMAKE_RANLIB ${TOOLCHAIN_PREFIX}-ranlib)

# target environment on the build host system
set (CMAKE_FIND_ROOT_PATH /usr/bin)

# modify default behavior of FIND_XXX() commands
set (CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set (CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set (CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# Custom/own flags only understood by this project's CMakeLists.txt
set (
    CMAKE_CXX_LINK_FLAGS_TOOLCHAIN
    "-static-libstdc++ -static-libgcc -Wl,-Bstatic, -lwinpthread -Wl,-Bdynamic"
    )

# workarounds
set (
    CMAKE_CXX_FLAGS_TOOLCHAIN
    -DJUCE_USE_WINDOWS_MEDIA_FORMAT=0
    )
