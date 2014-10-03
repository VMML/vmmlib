# Copyright (c) BBP/EPFL 2011-2014, Stefan.Eilemann@epfl.ch
# Change this number when adding tests to force a CMake run: 2

if(NOT Boost_FOUND)
  return()
endif()

set(TEST_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
include(CommonCTest)

include(InstallFiles)
install_files(share/vmmlib/tests FILES ${TEST_FILES} COMPONENT examples)
