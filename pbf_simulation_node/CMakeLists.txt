cmake_minimum_required(VERSION 3.26)
project(PBFNode)

set(CMAKE_CXX_STANDARD 17)

set(TARGET PBFNode)

# Load Source
add_library(
        ${TARGET}
        SHARED
        Entrance.cpp
        solver.cpp
        solver.h
        SpatialGrid.h
        SpatialGrid.hpp
        sph_kernel.cpp
        sph_kernel.h
)

# Link Houdini Toolkit
target_link_libraries(
        ${TARGET}
        PUBLIC
        Houdini
)
target_include_directories(
        ${TARGET}
        PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}
)
houdini_configure_target(${TARGET})
