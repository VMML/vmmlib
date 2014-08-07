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

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/CMake/common)
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

if(LAPACK_FOUND)
  add_definitions(-DVMMLIB_USE_LAPACK)
  set(LAPACK_HEADERS
    vmmlib/lapack/detail/clapack.h
    vmmlib/lapack/detail/f2c.h )
endif()

if(CBLAS_FOUND)
  add_definitions(-DVMMLIB_USE_BLAS)
endif()

add_subdirectory(tests)

install(FILES ${HEADERS} DESTINATION include/vmmlib COMPONENT dev)
install(FILES ${DOCS} DESTINATION share/vmmlib COMPONENT dev)

if(LAPACK_FOUND)
  install(FILES ${LAPACK_HEADERS} DESTINATION include/vmmlib/lapack/detail
    COMPONENT dev)
endif()

include(DoxygenRule) # must be after all targets
include(CPackConfig)
include(CTest)
include(PackageConfig)
