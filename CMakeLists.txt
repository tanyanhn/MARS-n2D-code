cmake_minimum_required(VERSION 3.22)
# clang is faster than gcc, O2 is faster than O3
include(${MARS_COMPILER_FILE})
# message("MARS_COMPILER_FILE: ${MARS_COMPILER_FILE}")
project(MARS_n2D)

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES "${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES}")
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT DIM)
  set(DIM "2")
endif()

# Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -DDIM=${DIM}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")
set(CMAKE_CXX_STANDARD 20)

message("System : ${CMAKE_SYSTEM_NAME}")
message("Build type : ${CMAKE_BUILD_TYPE}")
message("CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")


# LINK LIBS
find_package(PkgConfig REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(GSL REQUIRED)  
find_package(Eigen3 CONFIG REQUIRED)
find_package(CppUnit CONFIG REQUIRED)
pkg_check_modules(LAPACKE REQUIRED lapacke)
# message("lib: " ${Eigen_LIBRARY})

option(USE_OPENMP "Whether to use the OpenMP parallelism." ON)
if (USE_OPENMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENMP -fopenmp")
endif()

# INCLUDE 
include_directories(src)

# MAIN LIB
set(MAIN_LIB "MARS")
add_subdirectory(src)


# TEST
option(ENABLE_TEST "Whether to enable the tests. " ON)
if (ENABLE_TEST)
  enable_testing()
  add_subdirectory(test)
endif()


# FORMAT
file(GLOB_RECURSE FORMAT_FILES CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/*.[ch]pp
    ${CMAKE_SOURCE_DIR}/*.[ch]
)
add_custom_target(format
COMMAND clang-format -i --style=file ${FORMAT_FILES}
COMMENT "Formatting source files with clang-format"
VERBATIM
)
