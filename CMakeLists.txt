cmake_minimum_required(VERSION 3.5)
project(bays)

if(NOT DIM)
    set(DIM "2")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -DDIM=${DIM}")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_STANDARD 17)

option(USE_OPENMP "Whether to use the OpenMP parallelism." OFF)
if (USE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENMP -fopenmp")
endif()

message(STATUS "System : ${CMAKE_SYSTEM_NAME}")
message(STATUS "Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE : ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG : ${CMAKE_CXX_FLAGS_DEBUG}")

find_package(PkgConfig REQUIRED)
pkg_check_modules(BLAS REQUIRED blas)
pkg_check_modules(CBLAS cblas)
pkg_check_modules(LAPACK REQUIRED lapack)
pkg_check_modules(LAPACKE REQUIRED lapacke)
pkg_check_modules(CPPUNIT REQUIRED cppunit)
include_directories(${BLAS_INCLUDE_DIRS} ${CBLAS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS} ${LAPACKE_INCLUDE_DIRS} ${CPPUNIT_INCLUDE_DIRS})
link_directories(${BLAS_LIBRARY_DIRS} ${CBLAS_LIBRARY_DIRS} ${LAPACK_LIBRARY_DIRS} ${LAPACKE_LIBRARY_DIRS} ${CPPUNIT_LIBRARY_DIRS})

enable_testing()
include_directories(src)
add_subdirectory(src)
add_subdirectory(test)



