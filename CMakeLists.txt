cmake_minimum_required(VERSION 3.16)

project(
    prand
    VERSION 0.0.1
    DESCRIPTION "A C99 parallel Mersenne Twister implementation"
    HOMEPAGE_URL "https://github.com/cheng-zhao/prand/"
    LANGUAGES C
)

# C99 required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# output all artifacts in top-level build directory. multi-config generators
# will have an extra subdirectory per config
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# determine if generator is multi-config or not
get_property(PRAND_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
# use Debug config by default. can set for single-config generators with
# -DCMAKE_BUILD_TYPE=Release on the command line
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# build examples
option(PRAND_BUILD_EXAMPLES "Build example programs" ON)
# if building static, build using position-independent code
option(PRAND_STATIC_PIC "BUild static library as PIC" OFF)

# print config options
if(PRAND_IS_MULTI_CONFIG)
    # multi-config generator allows config selection at build time
    message(STATUS "Build config: Multi")
else()
    message(STATUS "Build config: ${CMAKE_BUILD_TYPE}")
endif()
if(BUILD_SHARED_LIBS)
    message(STATUS "Build libraries: Shared")
else()
    message(STATUS "Build libraries: Static")
endif()
if(PRAND_BUILD_EXAMPLES)
    message(STATUS "Build examples: Yes")
else()
    message(STATUS "Build examples: No")
endif()
if(PRAND_STATIC_PIC)
    message(STATUS "Build PIC static: Yes")
else()
    message(STATUS "Build PIC static: No")
endif()

# compile options
if(MSVC)
    add_compile_options(
        /Wall
        # C4710: function not inlined
        /wd4710
        # C4711: function selected for automatic inline
        /wd4711
        # C5045: compiler will insert Spectre mitigation with /QSpectre
        /wd5045
    )
else()
    add_compile_options(-Wall)
endif()

# enable add_test command
include(CTest)

# examples
if(PRAND_BUILD_EXAMPLES)
    add_subdirectory(example)
endif()

add_subdirectory(src)
