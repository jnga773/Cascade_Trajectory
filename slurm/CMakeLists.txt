# executable to use (see if better way than specifying here)
set(ARRAY_JOB_EXECUTABLE "${CMAKE_BINARY_DIR}/BMtest/cctest_matrix")

message(STATUS "Preparing array job")
message(STATUS "Array job executable: ${ARRAY_JOB_EXECUTABLE}")

# copy files in build directory
configure_file(run-array.py ${CMAKE_BINARY_DIR} @ONLY)
configure_file(_run-array.sl ${CMAKE_BINARY_DIR} @ONLY)
