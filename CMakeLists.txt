# Project name
project (PoreDy)

# CMake version minimum required
cmake_minimum_required (VERSION 2.6)

# Include Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

include_directories(./include/Eigen)

# Add source
aux_source_directory(./source DIR_SRCS)

add_executable(PoreDy ${DIR_SRCS})

target_link_libraries (PoreDy Eigen3::Eigen)