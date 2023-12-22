cmake_minimum_required(VERSION 3.24)
project(fluid_houdini_test)

set(CMAKE_CXX_STANDARD 14)

# Find Houdini
set(Houdini_PATH "E:/Program Files/Side Effects Software/Houdini 20.0.506")
set(Houdini_DIR ${Houdini_PATH}/toolkit/cmake)
find_package(Houdini REQUIRED)

add_subdirectory(basicSPH_simulation_node)
add_subdirectory(basicSPH_simulation_normal)
add_subdirectory(PCISPH_simulation_node)