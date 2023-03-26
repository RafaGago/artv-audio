add_library (fonts STATIC)
target_include_directories (fonts PUBLIC ${PROJECT_SOURCE_DIR}/resources)
# each font is a translation unit; the linker should delete the unreferenced ones
target_sources(
  fonts
  PRIVATE
  OpenSans-Bold.ttf.cpp
  PassionsConflict-Regular.ttf.cpp
  )
