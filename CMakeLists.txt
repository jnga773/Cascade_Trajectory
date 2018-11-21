cmake_minimum_required (VERSION 3.0)

enable_language(Fortran)

option(DOUBLE "Set precision to double" OFF)

if (DOUBLE)
	set(DOUBLE_OPT "-r8")
	if (CMAKE_COMPILE_ID STREQUAL "GNU")
		set(DOUBLE_OPT "-fdefault-real-8")
	endif()
endif()

add_executable(cross_correlation cross_correlation.f90)

