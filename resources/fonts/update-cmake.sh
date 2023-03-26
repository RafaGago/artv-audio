#!/bin/bash

SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)
cd "$SCRIPT_DIR"

rm -f CMakeLists.txt
cat << EOF >> CMakeLists.txt
add_library (fonts STATIC)
target_include_directories (fonts PUBLIC \${PROJECT_SOURCE_DIR}/resources)
EOF

SOURCES=""
for CPP in $(ls *.cpp); do
  SOURCES="$SOURCES\n  $CPP"
done

echo "# each font is a translation unit; the linker should delete the unreferenced ones" >> CMakeLists.txt
printf "target_sources(\n  fonts\n  PRIVATE${SOURCES}\n  )\n" >> CMakeLists.txt
