# Project name ------------------------------------------------------------------------------------

project (PoreDy)

# Version number ----------------------------------------------------------------------------------

set (PoreDy_VERSION_MAJOR 0)

set (PoreDy_VERSION_MINOR 6)

# CMake version minimum required ------------------------------------------------------------------

cmake_minimum_required (VERSION 2.6)

# Include Eigen -----------------------------------------------------------------------------------

find_package(Eigen3 REQUIRED NO_MODULE)

include_directories(./include/Eigen)

# Add source --------------------------------------------------------------------------------------

aux_source_directory(./source DIR_SRCS)

add_executable(PoreDy ${DIR_SRCS})

target_link_libraries (PoreDy Eigen3::Eigen)

# Set compiling environment -----------------------------------------------------------------------

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -std=c++11 -O0 -Wall -g -ggdb")

set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -std=c++11 -O3 -Wall")