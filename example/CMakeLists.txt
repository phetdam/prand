cmake_minimum_required(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})

# multistream
add_executable(multistream multistream.c)
target_link_libraries(multistream PRIVATE prand)
# allow use as a test
add_test(NAME multistream COMMAND multistream)
