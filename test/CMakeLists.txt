add_library(Catch ${CMAKE_SOURCE_DIR}/test/catch_amalgamated.cpp)
# SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} $ENV{CXXFLAGS} -O0 -Wall
# -g2 -ggdb -fsanitize=address") SET(CMAKE_CXX_FLAGS_RELEASE
# "${CMAKE_CXX_FLAGS_RELEASE} $ENV{CXXFLAGS} -O2 -g")

if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(COMPLIE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  set(COMPLIE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

set(LIST_TEST testPolynomial TestCutCell TestInterfaceGraph TestMARSn2D)
foreach(test ${LIST_TEST})
  message(STATUS "test ${test}")
  add_executable(${test} ${test}.cpp)
  # Add the -std=c++11 flag as an example
  set_source_files_properties(${test}.cpp PROPERTIES COMPILE_FLAGS
                                                     "${COMPLIE_FLAGS}")
  # target_link_libraries(${test} ${MAIN_LIB} ${COMMON_LIB})
  target_link_libraries(
    ${test}
    ${MAIN_LIB}
    Catch
    fmt::fmt
    Eigen3::Eigen
    ${LAPACKE_LIBRARIES})
  add_test(
    NAME ${test}
    COMMAND ${test}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
  message("-- Target compiler flags are: ${COMPLIE_FLAGS}")
endforeach(test)

# set(LIST_TEST TestCurveFactory TestBA TestIT TestFL TestIT)
set(LIST_TEST)
foreach(test ${LIST_TEST})
  message(STATUS "test ${test}")
  add_executable(${test} ${test}.cpp)
  # Add the -std=c++11 flag as an example
  set_source_files_properties(${test}.cpp PROPERTIES COMPILE_FLAGS
                                                     "${COMPLIE_FLAGS}")
  # target_link_libraries(${test} ${MAIN_LIB} ${COMMON_LIB})
  target_link_libraries(${test} ${MAIN_LIB} fmt::fmt Eigen3::Eigen
                        ${LAPACKE_LIBRARIES})
  add_test(
    NAME ${test}
    COMMAND ${test}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
  message("-- Target compiler flags are: ${COMPLIE_FLAGS}")
endforeach(test)

add_executable(
  UnitTest
  UnitTest.cpp
  TestSegmentsIntersector.H
  TestSegmentsIntersector.cpp
  TestLocater.h
  TestLocater.cpp
  TestYinSetsBooleanOps.cpp
  TestYinSetsBooleanOps.H
  TestOrientedJordanCurve.h
  TestOrientedJordanCurve.cpp)
target_link_libraries(UnitTest ${CppUnit_LIBRARY} ${MAIN_LIB} fmt::fmt)

add_test(
  NAME UnitTest
  COMMAND UnitTest
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
