cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++17")

project(pro)
add_executable(bee main.cpp)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)
set(GMSH_LIBRARIES gmsh)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Укажем библиотеки для линковки
# ... тут только VTK
# target_link_libraries(cubic2d ${VTK_LIBRARIES})
# target_link_libraries(cubic3d ${VTK_LIBRARIES})

# ... а тут ещё и gmsh
target_link_libraries(bee ${VTK_LIBRARIES})
target_link_libraries(bee ${GMSH_LIBRARIES})