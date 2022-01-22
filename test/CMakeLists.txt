add_executable(TestCurveFactory TestCurveFactory.cpp)
target_link_libraries(TestCurveFactory bays)

#add_executable(TestBA TestBA.cpp)
#target_link_libraries(TestBA bays)

#add_executable(TestIT TestIT.cpp)
#target_link_libraries(TestIT bays)

add_executable(UnitTest UnitTest.cpp
        #TestSegmentsIntersector.H
        #TestSegmentsIntersector.cpp
        #TestLocater.h
        #TestLocater.cpp
        #TestYinSetsBooleanOps.cpp
        #TestYinSetsBooleanOps.H
        TestOrientedJordanCurve.h 
        TestOrientedJordanCurve.cpp)
target_link_libraries(UnitTest ${CPPUNIT_LIBRARIES} bays)

add_test(NAME UnitTest COMMAND UnitTest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)

