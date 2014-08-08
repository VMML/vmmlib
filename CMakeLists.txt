# Copyright (c) 2011 Daniel Pfeifer <daniel@pfeifer-mail.de>

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(vmmlib)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 7)
set(VERSION_PATCH 0)
set(VERSION_ABI 2)
set(LAST_RELEASE 1.7.0)

set(GITTARGETS_RELEASE_BRANCH minor)

option(VMMLIB_USE_BLAS "Disable to turn blas detection off" ON)
mark_as_advanced(VMMLIB_USE_BLAS)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake/common)
include(GitExternal)
include(Common)

include(Buildyard)
if(BUILDYARD_STOP)
  return()
endif()

include(FindPackages)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow -Wno-sign-compare -Wno-extra -Wno-strict-aliasing")
endif()

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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vmmlib/version.in.hpp
  ${OUTPUT_INCLUDE_DIR}/vmmlib/version.hpp @ONLY)

set(HEADERS
  ${OUTPUT_INCLUDE_DIR}/vmmlib/version.hpp
  vmmlib/aabb.hpp
  vmmlib/blas_daxpy.hpp
  vmmlib/blas_dgemm.hpp
  vmmlib/blas_dot.hpp
  vmmlib/blas_includes.hpp
  vmmlib/blas_types.hpp
  vmmlib/cp3_tensor.hpp
  vmmlib/enable_if.hpp
  vmmlib/exception.hpp
  vmmlib/frustum.hpp
  vmmlib/frustum_culler.hpp
  vmmlib/intersection.hpp
  vmmlib/jacobi_solver.hpp
  vmmlib/lapack.hpp
  vmmlib/lapack/detail/clapack.h
  vmmlib/lapack/detail/f2c.h
  vmmlib/lapack_gaussian_elimination.hpp
  vmmlib/lapack_includes.hpp
  vmmlib/lapack_linear_least_squares.hpp
  vmmlib/lapack_svd.hpp
  vmmlib/lapack_sym_eigs.hpp
  vmmlib/lapack_types.hpp
  vmmlib/linear_least_squares.hpp
  vmmlib/lowpass_filter.hpp
  vmmlib/math.hpp
  vmmlib/matrix.hpp
  vmmlib/matrix_functors.hpp
  vmmlib/matrix_pseudoinverse.hpp
  vmmlib/matrix_traits.hpp
  vmmlib/qr_decomposition.hpp
  vmmlib/qtucker3_tensor.hpp
  vmmlib/quaternion.hpp
  vmmlib/svd.hpp
  vmmlib/t3_converter.hpp
  vmmlib/t3_hooi.hpp
  vmmlib/t3_hopm.hpp
  vmmlib/t3_hosvd.hpp
  vmmlib/t3_ihooi.hpp
  vmmlib/t3_ihopm.hpp
  vmmlib/t3_ttm.hpp
  vmmlib/t4_converter.hpp
  vmmlib/t4_hooi.hpp
  vmmlib/t4_hosvd.hpp
  vmmlib/t4_ttm.hpp
  vmmlib/tensor3.hpp
  vmmlib/tensor3_iterator.hpp
  vmmlib/tensor4.hpp
  vmmlib/tensor_mmapper.hpp
  vmmlib/tensor_stats.hpp
  vmmlib/tucker3_exporter.hpp
  vmmlib/tucker3_importer.hpp
  vmmlib/tucker3_tensor.hpp
  vmmlib/tucker4_tensor.hpp
  vmmlib/util.hpp
  vmmlib/validator.hpp
  vmmlib/vector.hpp
  vmmlib/vector_traits.hpp
  vmmlib/visibility.hpp
  vmmlib/vmmlib.hpp
  vmmlib/vmmlib_config.hpp
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
    vmmlib/lapack/detail/clapack.h
    vmmlib/lapack/detail/f2c.h )
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

list(APPEND CPPCHECK_EXTRA_ARGS --suppress=*:${CMAKE_SOURCE_DIR}/tests/*)

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
  add_custom_target(vmmlib_tests COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS vmmlib_test)
  add_custom_target(vmmlib-test DEPENDS vmmlib_tests) # Travis CI hack
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
