

set(VMMLIB_PACKAGE_VERSION 1.7.0)
set(VMMLIB_REPO_URL https://github.com/Eyescale/vmmlib.git)
set(VMMLIB_DEPENDS OPTIONAL OpenMP Boost)
set(VMMLIB_BOOST_COMPONENTS "unit_test_framework")

if(CI_BUILD_COMMIT)
  set(VMMLIB_REPO_TAG ${CI_BUILD_COMMIT})
else()
  set(VMMLIB_REPO_TAG master)
endif()
set(VMMLIB_FORCE_BUILD ON)
set(VMMLIB_SOURCE ${CMAKE_SOURCE_DIR})