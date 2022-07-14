cmake_minimum_required(VERSION 3.0)
project(vecmath_test)
set(CMAKE_CONFIGURATION_TYPES "Release")
file(GLOB SRC *.cpp)
add_executable(vecmath_test ${SRC})

if (MSVC)
    add_compile_options(/O2)
else()
    add_compile_options(-O2)
endif()

target_include_directories(vecmath_test PUBLIC .)
target_link_libraries(vecmath_test m)
