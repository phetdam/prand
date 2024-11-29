cmake_minimum_required(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})

# libprand
add_library(prand mrg32k3a.c mt19937.c mt19937_poly.c prand.c)
# public header
set_target_properties(prand PROPERTIES PUBLIC_HEADER header/prand.h)
# use "d" suffix on Windows to distinguish use of debug C runtime
if(MSVC)
    set_target_properties(prand PROPERTIES OUTPUT_NAME_DEBUG prandd)
endif()
# install rule
install(TARGETS prand)
