
if(MSVC)
  return()
endif()

set(TESTS
  frustum_test.cpp
  intersection_test.cpp
  jacobi_test.cpp
  matrix_test.cpp
  perf_test.cpp
  performance_test.cpp
  qr_decomposition_test.cpp
  quaternion_test.cpp
  svd_test.cpp
  timer.cpp
  unit_test.cpp
  unit_test_globals.cpp
  lowpass_filter_test.cpp
  vector_test.cpp
  util_test.cpp
  vmmlib_unit_tests_main.cpp )

set(TEST_LIBRARIES)

if(LAPACK_FOUND)
  list(APPEND TEST_LIBRARIES ${LAPACK_LIBRARIES})
  list(APPEND TESTS
    lapack_gaussian_elimination_test.cpp
    lapack_linear_least_squares_test.cpp
    lapack_svd_test.cpp
    lapack_sym_eigs_test.cpp
    svd_lapack_vs_old.cpp )
endif()
if(CBLAS_FOUND)
  list(APPEND TEST_LIBRARIES ${CBLAS_LIBRARIES})
  list(APPEND TESTS
    blas_daxpy_test.cpp
    blas_dgemm_test.cpp
    blas_dot_test.cpp )
  if(LAPACK_FOUND)
    list(APPEND TESTS
      tensor3_iterator_test.cpp
      tensor3_test.cpp
      tucker3_exporter_importer_test.cpp
      tucker3_tensor_test.cpp
      t4_hooi_test.cpp
      t4_converter_test.cpp
      tensor4_test.cpp
      t4_ttm_test.cpp
      t3_ttm_test.cpp
      t3_hosvd_test.cpp
      t3_ihooi_test.cpp
      t3_hopm_test.cpp
      t3_ihopm_test.cpp
      t3_hooi_test.cpp
      qtucker3_tensor_test.cpp
      matrix_pseudoinverse_test.cpp
      cp3_tensor_test.cpp )
  endif()
endif()

list(APPEND CPPCHECK_EXTRA_ARGS --suppress=*:${CMAKE_CURRENT_SOURCE_DIR}/*)

add_executable(vmmlib_test ${TESTS})
get_target_property(EXECUTABLE vmmlib_test LOCATION)
STRING(REGEX REPLACE "\\$\\(.*\\)" "\${CTEST_CONFIGURATION_TYPE}"
  EXECUTABLE "${EXECUTABLE}")
add_test(vmmlib_test ${EXECUTABLE})
target_link_libraries(vmmlib_test ${TEST_LIBRARIES})

# workaround: 'make test' does not build tests beforehand
add_custom_target(tests COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS vmmlib_test)
add_custom_target(vmmlib-test DEPENDS tests) # Travis CI hack
