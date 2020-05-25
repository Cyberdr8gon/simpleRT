cmake_minimum_required ( VERSION 3.12.0 )

project ( SimpleRT LANGUAGES CXX)

# Set to c++17
set ( CMAKE_CXX_STANDARD 17 )
list(APPEND CMAKE_MODULE_PATH "$ENV{Eigen3_ROOT}/cmake")

SET(Eigen3_DIR "$ENV{Eigen3_ROOT}" CACHE LOCATION "Where is Eigen3?")
find_package(Eigen3 3.3 REQUIRED)

include_directories("third_party/stb")

message(STATUS "Eigen Dir " ${EIGEN3_INCLUDE_DIR})

add_definitions(${EIGEN3_DEFINITIONS})
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories("include")

add_executable(simpleRT "src/main.cpp")
