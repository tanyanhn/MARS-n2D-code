add_subdirectory(Core)
add_subdirectory(YinSet)
add_subdirectory(InterfaceTracking)
add_subdirectory(GeometricFlow)
# add_subdirectory(Marsn2D)
add_subdirectory(Recorder)

add_library(${MAIN_LIB} STATIC ${Core_Files} ${YinSet_Files} ${InterfaceTracking_Files} ${GeometricFlow_Files}
${Recorder_Files} ${Marsn2D_Files})
target_link_libraries(${MAIN_LIB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES})
