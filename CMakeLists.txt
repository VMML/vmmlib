# Copyright (c) 2011 Daniel Pfeifer <daniel@pfeifer-mail.de>
#               2011-2014 Stefan.Eilemann@epfl.ch

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(vmmlib)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 8)
set(VERSION_PATCH 0)
set(VERSION_ABI 3)

configure_file(${PROJECT_SOURCE_DIR}/vmmlib/version.in.hpp
               ${PROJECT_SOURCE_DIR}/vmmlib/version.hpp @ONLY)

set(VMMLIB_COMPILE_TESTS TRUE CACHE BOOL "Compile 'tests' and 'tests_old' directories")

if(VMMLIB_COMPILE_TESTS)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
                                ${CMAKE_SOURCE_DIR}/CMake/common)
  include(Common)
  include(FindPackages)

  common_package(Boost COMPONENTS unit_test_framework)

  add_subdirectory(tests)
  #add_subdirectory(tests_old)
endif()
