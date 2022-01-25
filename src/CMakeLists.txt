add_subdirectory(Core)
add_subdirectory(YinSet)
# add_subdirectory(InterfaceTracking)
# add_subdirectory(GeometricFlow)

add_library(bays STATIC ${Core_Files} ${YinSet_Files} ${InterfaceTracking_Files} ${GeometricFlow_Files})
target_link_libraries(bays ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES} gfortran)
