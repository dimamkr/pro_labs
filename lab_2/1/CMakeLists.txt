cmake_minimum_required(VERSION 3.21)

set(PROJECT_NAME demo_poisson)
project(${PROJECT_NAME} LANGUAGES C CXX)

# Если мы в conda окружении, используем его как корень для поиска
if(DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_FIND_ROOT_PATH "$ENV{CONDA_PREFIX}")
    set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
    set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
    set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
    set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

    # Указываем HDF5 явно
    set(HDF5_ROOT "$ENV{CONDA_PREFIX}")
    set(HDF5_PREFER_PARALLEL TRUE)
    set(HDF5_DIR "$ENV{CONDA_PREFIX}/lib/cmake/hdf5")
endif()

if(NOT TARGET dolfinx)
  find_package(DOLFINX REQUIRED)
endif()

include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
check_symbol_exists(PETSC_USE_COMPLEX petscsystypes.h PETSC_SCALAR_COMPLEX)
check_symbol_exists(PETSC_USE_REAL_DOUBLE petscsystypes.h PETSC_REAL_DOUBLE)

# Add target to compile UFL files
if(PETSC_SCALAR_COMPLEX EQUAL 1)
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=complex128")
  else()
    set(SCALAR_TYPE "--scalar_type=complex64")
  endif()
else()
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=float64")
  else()
    set(SCALAR_TYPE "--scalar_type=float32")
  endif()
endif()
add_custom_command(
  OUTPUT poisson.c
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/poisson.py ${SCALAR_TYPE}
  VERBATIM
  DEPENDS poisson.py
  COMMENT "Compile poisson.py using FFCx"
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(${PROJECT_NAME} main.cpp ${CMAKE_CURRENT_BINARY_DIR}/poisson.c)
target_link_libraries(${PROJECT_NAME} dolfinx)

# Set C++20 standard
set(CMAKE_CXX_EXTENSIONS OFF)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

# Убираем -Werror, чтобы не ломать сборку из-за предупреждений
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-comment" HAVE_NO_MULTLINE)
if(HAVE_NO_MULTLINE)
  target_compile_options(${PROJECT_NAME} PRIVATE -Wno-comment)
endif()
# Можно добавить -Wall -Wextra, но без -Werror
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)

# Test targets
set(TEST_PARAMETERS2 -np 2 "./${PROJECT_NAME}")
set(TEST_PARAMETERS3 -np 3 "./${PROJECT_NAME}")
add_test(NAME ${PROJECT_NAME}_mpi_2 COMMAND "mpirun" ${TEST_PARAMETERS2})
add_test(NAME ${PROJECT_NAME}_mpi_3 COMMAND "mpirun" ${TEST_PARAMETERS3})
add_test(NAME ${PROJECT_NAME}_serial COMMAND ${PROJECT_NAME})