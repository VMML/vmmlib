# Copyright (c) 2012 Stefan Eilemann <eile@eyescale.ch>

# Info: http://www.itk.org/Wiki/CMake:Component_Install_With_CPack

set(CPACK_PACKAGE_VENDOR "vmmlib.sourceforge.net")
set(CPACK_PACKAGE_CONTACT "Stefan Eilemann <eile@eyescale.ch>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "Templatized C++ vector and matrix math library")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6")

set(CPACK_MACPORTS_CATEGORY devel)

if(RELEASE_VERSION)
  set(CPACK_COMPONENTS_ALL dev)
else()
  set(CPACK_COMPONENTS_ALL unspecified dev)
endif()

set(CPACK_COMPONENT_UNSPECIFIED_DISPLAY_NAME "Unspecified")
set(CPACK_COMPONENT_UNSPECIFIED_DESCRIPTION
  "Unspecified Component - set COMPONENT in CMake install() command")

set(CPACK_COMPONENT_DEV_DISPLAY_NAME
  "${CPACK_PROJECT_NAME} Development Files")
set(CPACK_COMPONENT_DEV_DESCRIPTION
  "Header and Library Files for ${CPACK_PROJECT_NAME} Development")

include(CommonCPack)
