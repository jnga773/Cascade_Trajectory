cmake_minimum_required (VERSION 3.0)
project(cross_correlation LANGUAGES Fortran)

# add testing, can be disabled with -DBUILD_TESTING=OFF
include(CTest)

# default to a release (optimised) build
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# make sure the Cray compiler uses -O3 in RELEASE
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
endif()

# promotes real to double precision
option(DOUBLE_PRECISION "Set precision to double" ON)
if (DOUBLE_PRECISION)
	set(DOUBLE_OPT "-r8")
	if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		set(DOUBLE_OPT "-fdefault-real-8")
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        set(DOUBLE_OPT "-s real64")
	endif()
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DOUBLE_OPT}")
endif()

add_executable(cross_correlation cross_correlation.f90)

# add BMtest directory
add_subdirectory(BMtest)

# add slurm directory
add_subdirectory(slurm)
