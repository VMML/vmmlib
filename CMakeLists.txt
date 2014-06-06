# Copyright (c) 2011 Daniel Pfeifer <daniel@pfeifer-mail.de>

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(vmmlib)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 7)
set(VERSION_PATCH 0)
set(RELEASE_VERSION OFF) # OFF or ABI version
set(LAST_RELEASE 1.6.1)

option(VMMLIB_USE_BLAS "Disable to turn blas detection off" ON)
mark_as_advanced(VMMLIB_USE_BLAS)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/CMake/common)
include(GitExternal)
include(Common)

include(Buildyard)
if(BUILDYARD_STOP)
  return()
endif()

set(MACPORTS_DIR eyescalePorts) # override for now...

include(FindPackages)
include(CheckLibraryExists)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow -Wno-sign-compare -Wno-extra -Wno-strict-aliasing")
endif()

include_directories(BEFORE ${CMAKE_SOURCE_DIR}/include)

if(VMMLIB_USE_BLAS)
    if(APPLE)
      find_library(LAPACK_LIBRARIES Accelerate)
      if(LAPACK_LIBRARIES)
        set(LAPACK_LIBRARIES lapack)
        set(CBLAS_LIBRARIES blas)
      endif(LAPACK_LIBRARIES)
    endif()
endif()

set(DOCS README.md AUTHORS LICENSE.txt ACKNOWLEDGEMENTS)

update_file(${CMAKE_CURRENT_SOURCE_DIR}/include/vmmlib/version.in.hpp
  ${OUTPUT_INCLUDE_DIR}/vmmlib/version.hpp)

set(HEADERS
  ${OUTPUT_INCLUDE_DIR}/vmmlib/version.hpp
  include/vmmlib/blas_daxpy.hpp
  include/vmmlib/blas_dgemm.hpp
  include/vmmlib/blas_dot.hpp
  include/vmmlib/blas_includes.hpp
  include/vmmlib/blas_types.hpp
  include/vmmlib/aabb.hpp
  include/vmmlib/cp3_tensor.hpp
  include/vmmlib/enable_if.hpp
  include/vmmlib/exception.hpp
  include/vmmlib/frustum.hpp
  include/vmmlib/frustum_culler.hpp
  include/vmmlib/intersection.hpp
  include/vmmlib/jacobi_solver.hpp
  include/vmmlib/lapack.hpp
  include/vmmlib/lapack/detail/clapack.h
  include/vmmlib/lapack/detail/f2c.h
  include/vmmlib/lapack_gaussian_elimination.hpp
  include/vmmlib/lapack_includes.hpp
  include/vmmlib/lapack_linear_least_squares.hpp
  include/vmmlib/lapack_svd.hpp
  include/vmmlib/lapack_sym_eigs.hpp
  include/vmmlib/lapack_types.hpp
  include/vmmlib/linear_least_squares.hpp
  include/vmmlib/lowpass_filter.hpp
  include/vmmlib/math.hpp
  include/vmmlib/matrix.hpp
  include/vmmlib/matrix_functors.hpp
  include/vmmlib/matrix_pseudoinverse.hpp
  include/vmmlib/matrix_traits.hpp
  include/vmmlib/qr_decomposition.hpp
  include/vmmlib/qtucker3_tensor.hpp
  include/vmmlib/quaternion.hpp
  include/vmmlib/svd.hpp
  include/vmmlib/t3_converter.hpp
  include/vmmlib/t3_hooi.hpp
  include/vmmlib/t3_hopm.hpp
  include/vmmlib/t3_hosvd.hpp
  include/vmmlib/t3_ihopm.hpp
  include/vmmlib/t3_ihooi.hpp
  include/vmmlib/t3_ttm.hpp
  include/vmmlib/tensor_mmapper.hpp
  include/vmmlib/tensor_stats.hpp
  include/vmmlib/tensor3.hpp
  include/vmmlib/tensor3_iterator.hpp
  include/vmmlib/tucker3_exporter.hpp
  include/vmmlib/tucker3_importer.hpp
  include/vmmlib/tucker3_tensor.hpp
  include/vmmlib/util.hpp
  include/vmmlib/validator.hpp
  include/vmmlib/vector.hpp
  include/vmmlib/vector_traits.hpp
  include/vmmlib/visibility.hpp
  include/vmmlib/vmmlib.hpp
  include/vmmlib/vmmlib_config.hpp
  include/vmmlib/tensor4.hpp
  include/vmmlib/t4_converter.hpp
  include/vmmlib/t4_hooi.hpp
  include/vmmlib/t4_hosvd.hpp
  include/vmmlib/t4_ttm.hpp
  include/vmmlib/tucker4_tensor.hpp
)

set(TESTS
  tests/frustum_test.cpp
  tests/intersection_test.cpp
  tests/jacobi_test.cpp
  tests/matrix_test.cpp
  tests/perf_test.cpp
  tests/performance_test.cpp
  tests/qr_decomposition_test.cpp
  tests/quaternion_test.cpp
  tests/svd_test.cpp
  tests/timer.cpp
  tests/unit_test.cpp
  tests/unit_test_globals.cpp
  tests/lowpass_filter_test.cpp
  tests/vector_test.cpp
  tests/util_test.cpp
  tests/vmmlib_unit_tests_main.cpp )

if(LAPACK_FOUND)
  add_definitions(-DVMMLIB_USE_LAPACK)
  set(LAPACK_HEADERS
    include/vmmlib/lapack/detail/clapack.h
    include/vmmlib/lapack/detail/f2c.h )
  list(APPEND TESTS
    tests/lapack_gaussian_elimination_test.cpp
    tests/lapack_linear_least_squares_test.cpp
    tests/lapack_svd_test.cpp
    tests/lapack_sym_eigs_test.cpp
    tests/svd_lapack_vs_old.cpp )
endif()

if(CBLAS_FOUND)
  add_definitions(-DVMMLIB_USE_BLAS)
  list(APPEND TESTS
    tests/blas_daxpy_test.cpp
    tests/blas_dgemm_test.cpp
    tests/blas_dot_test.cpp )
endif()

if(CBLAS_FOUND AND LAPACK_FOUND)
  add_definitions(-DVMMLIB_USE_BLAS)
  list(APPEND TESTS
    tests/tensor3_iterator_test.cpp
    tests/tensor3_test.cpp
    tests/tucker3_exporter_importer_test.cpp
    tests/tucker3_tensor_test.cpp
    tests/t4_hooi_test.cpp
    tests/t4_converter_test.cpp
    tests/tensor4_test.cpp
    tests/t4_ttm_test.cpp
    tests/t3_ttm_test.cpp
    tests/t3_hosvd_test.cpp
    tests/t3_ihooi_test.cpp
    tests/t3_hopm_test.cpp
    tests/t3_ihopm_test.cpp
    tests/t3_hooi_test.cpp
    tests/qtucker3_tensor_test.cpp
    tests/matrix_pseudoinverse_test.cpp
    tests/cp3_tensor_test.cpp )
endif()

list(APPEND CPPCHECK_EXTRA_ARGS
  --suppress=*:${CMAKE_SOURCE_DIR}/tests/*)

if(NOT MSVC)
  add_executable(vmmlib_test ${HEADERS} ${TESTS})
  if(LAPACK_FOUND)
    target_link_libraries(vmmlib_test ${LAPACK_LIBRARIES})
  endif()
  if(CBLAS_FOUND)
    target_link_libraries(vmmlib_test ${CBLAS_LIBRARIES})
  endif()

  get_target_property(EXECUTABLE vmmlib_test LOCATION)
  STRING(REGEX REPLACE "\\$\\(.*\\)" "\${CTEST_CONFIGURATION_TYPE}"
    EXECUTABLE "${EXECUTABLE}")
  add_test(vmmlib_test ${EXECUTABLE})

  # workaround: 'make test' does not build tests beforehand
  add_custom_target(tests COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS vmmlib_test)
  add_custom_target(vmmlib-test DEPENDS tests) # Travis CI hack
endif()

install(FILES ${HEADERS} DESTINATION include/vmmlib COMPONENT dev)
install(FILES ${DOCS} DESTINATION share/vmmlib COMPONENT dev)

if(LAPACK_FOUND)
   install(FILES ${LAPACK_HEADERS} DESTINATION include/vmmlib/lapack/detail COMPONENT dev)
endif()

include(DoxygenRule) # must be after all targets
include(CPackConfig)
include(CTest)
include(PackageConfig)
