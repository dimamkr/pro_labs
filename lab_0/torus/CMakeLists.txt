cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(gmshdemo)

add_executable(main main.cpp)

# --- Подключаем gmsh
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)
set(GMSH_LIBRARIES gmsh)

# --- Закончили с gmsh

# Укажем библиотеки для линковки
target_link_libraries(main ${GMSH_LIBRARIES})
