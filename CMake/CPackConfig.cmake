# Copyright (c) 2012 Stefan Eilemann <eile@eyescale.ch>

# Info: http://www.itk.org/Wiki/CMake:Component_Install_With_CPack

set(CPACK_PACKAGE_VENDOR "vmmlib.sourceforge.net")
set(CPACK_PACKAGE_CONTACT "Stefan Eilemann <eile@eyescale.ch>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "Templatized C++ vector and matrix math library")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
  "libstdc++6")

include(CommonCPack)
