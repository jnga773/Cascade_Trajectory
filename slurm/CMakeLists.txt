message(STATUS "Preparing array job")
message(STATUS "Array job executable: ${ARRAY_JOB_EXE}")

# copy files in build directory
configure_file(run-array.py ${CMAKE_BINARY_DIR} @ONLY)
configure_file(_run-array.sl ${CMAKE_BINARY_DIR} @ONLY)
