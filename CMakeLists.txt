cmake_minimum_required (VERSION 3.0)

enable_language(Fortran)

option(DOUBLE_PRECISION "Set precision to double" OFF)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if (DOUBLE_PRECISION)
	# promotes real to double precision
	set(DOUBLE_OPT "-r8")
	if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		set(DOUBLE_OPT "-fdefault-real-8")
	endif()
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DOUBLE_OPT}")
endif()

add_executable(cross_correlation cross_correlation.f90)

