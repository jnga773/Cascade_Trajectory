#
# CMake script to run a single test case
#
message(STATUS "Running trajectories test")
message(STATUS "  Test run directory: ${TEST_RUN_DIR}")
message(STATUS "  Test src directory: ${TEST_SRC_DIR}")
message(STATUS "  Test binary: ${TEST_BINARY}")
message(STATUS "  Test steps: ${TEST_STEPS}")
message(STATUS "  Compare script: ${COMPARE_SCRIPT}")
message(STATUS "  Reference file: ${REF_FILE}")
message(STATUS "  Python found: ${PYTHON_FOUND}")
message(STATUS "  Python: ${PYTHON}")

#
# make the test directory
#
execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${TEST_RUN_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${TEST_RUN_DIR})

#
# run the code
#
execute_process(
    COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR} ${TEST_BINARY} ${TEST_STEPS}
    RESULT_VARIABLE status
)
if (status)
    message(FATAL_ERROR "Error running BMtest code: '${status}'")
endif (status)

#
# check the results
#

# fail if no python is available
if (NOT PYTHON_FOUND)
    message(FATAL_ERROR "Failing test due to no Python available for comparing results")
endif()

# also fail if no reference file is available
if (NOT EXISTS ${REF_FILE})
    message(FATAL_ERROR "Failing test due to no reference file for comparing results: ${REF_FILE}")
endif()

# run the comparison script
execute_process(
    COMMAND ${CMAKE_COMMAND} -E chdir ${TEST_RUN_DIR}
        ${PYTHON} ${COMPARE_SCRIPT} ${REF_FILE}
    RESULT_VARIABLE status
)
if (status)
    message(FATAL_ERROR "Output files do not match: '${status}'")
endif (status)
