
if(MSVC)
  return()
endif()

set(TESTS
  intersection_test.cpp
  matrix_test.cpp
  perf_test.cpp
  performance_test.cpp
  qr_decomposition_test.cpp
  quaternion_test.cpp
  timer.cpp
  unit_test.cpp
  unit_test_globals.cpp
  lowpass_filter_test.cpp
  vector_test.cpp
  util_test.cpp
  vmmlib_unit_tests_main.cpp)

set(TEST_LIBRARIES)

list(APPEND CPPCHECK_EXTRA_ARGS --suppress=*:${PROJECT_SOURCE_DIR}/*)

add_executable(vmmlib_test ${TESTS})
if(CMAKE_VERSION VERSION_GREATER 2.8.10)
  add_test(NAME vmmlib_test COMMAND $<TARGET_FILE:vmmlib_test>)
else()
  get_target_property(EXECUTABLE vmmlib_test LOCATION)
  string(REGEX REPLACE "\\$\\(.*\\)" "\${CTEST_CONFIGURATION_TYPE}"
               EXECUTABLE "${EXECUTABLE}")
  add_test(vmmlib_test ${EXECUTABLE})
endif()
target_link_libraries(vmmlib_test ${TEST_LIBRARIES})

# workaround: 'make test' does not build tests beforehand
add_custom_target(old_tests COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS vmmlib_test)
