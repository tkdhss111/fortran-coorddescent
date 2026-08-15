cmake_minimum_required ( VERSION 3.16 )
project ( fortran-coorddescent LANGUAGES Fortran )

set ( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod )
if ( NOT CMAKE_BUILD_TYPE )
  set ( CMAKE_BUILD_TYPE Release )
endif ()

#
# Coarrays. -coarray=shared launches images through MPI Hydra, so the Intel
# oneAPI setvars.sh must be sourced before RUNNING, not only before building.
#
# Profile with -coarray=single: a shared build is already running N images, so
# a profile of it measures the wrong thing.
#
set ( COARRAY_MODE "shared" CACHE STRING "single | shared | distributed" )
set ( COARRAY_IMAGES "4" CACHE STRING "images for -coarray=shared" )

if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  set ( CA_FLAGS "-coarray=${COARRAY_MODE}" )
  if ( COARRAY_MODE STREQUAL "shared" )
    set ( CA_FLAGS "${CA_FLAGS} -coarray-num-images=${COARRAY_IMAGES}" )
  endif ()
else ()
  message ( FATAL_ERROR "Coarray support here is Intel-specific; build with ifx." )
endif ()

set ( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
set ( CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -check all -traceback" )
set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CA_FLAGS}" )

add_library ( coorddescent STATIC src/coord_descent_mo.f90 )

enable_testing ()
add_executable ( test_coord_descent test/test_coord_descent.f90 )
target_link_libraries ( test_coord_descent PRIVATE coorddescent )
add_test ( NAME test_coord_descent COMMAND test_coord_descent )
