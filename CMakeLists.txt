cmake_minimum_required(VERSION 3.16)

project(multifebe VERSION 2.0.1) # The version number must be changed here for every new change
                                 # Note that it must be changed manually too in src/version.f90
                                 # This is something that could be automated in the future

enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-std=gnu -ffree-line-length-none -Wall -Wextra -cpp -fimplicit-none -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -fopenmp -march=core2 -O3 -flto")
    set(bounds "-std=gnu -ffree-line-length-none -Wall -Wextra -cpp -fimplicit-none -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -fno-inline -Og")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

file(GLOB_RECURSE sources  src/*.f90 lib/fbem/src/*.f90)

add_executable(multifebe ${sources})

include_directories(lib/fbem/src/resources_shape_functions)
include_directories(lib/fbem/src/resources_quad_rules)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)


#
# Build and link to SDF
#
# add_custom_target(libsdf.a COMMAND make
                   # WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/sdf-0.75-RC4
# )
# target_link_libraries(multifebe ${CMAKE_CURRENT_SOURCE_DIR}/lib/sdf-0.75-RC4/libsdf.a)
# add_dependencies(multifebe libsdf.a)


#
# Link to OpenBLAS
#
set(BLA_VENDOR OpenBLAS)
set(BLA_STATIC ON)
find_package(BLAS REQUIRED)
if (UNIX)
  target_link_libraries(multifebe ${BLAS_LIBRARIES})
endif (UNIX)
#
# In Windows MSYS2 MINGW X64 installation, static linking of OpenBLAS fails.
# Despite BLA_STATIC is ON, it selects the dynamic library. We have
# to delete the ".dll" part of the filename in order to take the static one.
#
if (WIN32)
  string(REPLACE ".dll" "" BLAS_LIBRARIES_MOD ${BLAS_LIBRARIES})  
  message(STATUS "Modification for taking OpenBLAS static library ${BLAS_LIBRARIES_MOD}")
  target_link_libraries(multifebe ${BLAS_LIBRARIES_MOD})
endif (WIN32)

#
# For installing the program and documentation:
# $ make install
include(GNUInstallDirs)
install(TARGETS multifebe RUNTIME COMPONENT applications)
install(FILES "${CMAKE_SOURCE_DIR}/docs/manual.pdf" TYPE DOC COMPONENT documentation)
install(FILES "${CMAKE_SOURCE_DIR}/docs/gid/multifebe_unix.bas" TYPE DOC COMPONENT documentation)
install(FILES "${CMAKE_SOURCE_DIR}/docs/gid/multifebe_win.bas" TYPE DOC COMPONENT documentation)

#
# For generating an installer on GNU/Linux (debian)
# $ cpack -G "DEB" -C CPackConfig.cmake ..
#
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_CONTACT "jacobdavid.rodriguezbordon@ulpgc.es")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Multi-domain Finite Elements - Boundary Elements linear mechanics solver")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/mmc-siani-es/MultiFEBE")
set(CPACK_PACKAGE_VENDOR "Universidad de Las Palmas de Gran Canaria - SIANI-MMCE")
set(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/docs/img\\\\logo.png")
set(CPACK_COMPONENTS_ALL applications documentation)
set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "MultiFEBE")
set(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
if (UNIX)
  set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://github.com/mmc-siani-es/MultiFEBE")
  set(CPACK_GENERATOR "DEB;STGZ")
endif (UNIX)



#
# For generating an installer on Windows (in MSYS2 MINGW64 terminal)
# cpack -G "NSIS64" ..
# 
if (WIN32)
  set(CPACK_NSIS_DISPLAY_NAME "MultiFEBE")
  set(CPACK_NSIS_PACKAGE_NAME "MultiFEBE")
  set(CPACK_NSIS_CONTACT "jacobdavid.rodriguezbordon@ulpgc.es")
  set(CPACK_NSIS_URL_INFO_ABOUT "https://github.com/mmc-siani-es/MultiFEBE")
  set(CPACK_NSIS_HELP_LINK "https://github.com/mmc-siani-es/MultiFEBE")
  set(CPACK_NSIS_WELCOME_TITLE "Welcome to the MultiFEBE installer")
  # Installer welcome/unwelcome bitmap must be 164 x 314 px
  set(CPACK_NSIS_MUI_WELCOMEFINISHPAGE_BITMAP "${CMAKE_CURRENT_SOURCE_DIR}/docs/img\\\\installer_welcome.bmp")
  set(CPACK_NSIS_MUI_UNWELCOMEFINISHPAGE_BITMAP "${CMAKE_CURRENT_SOURCE_DIR}/docs/img\\\\installer_welcome.bmp")
  # Installer header bitmap must be 150 x 57 px
  set(CPACK_NSIS_MUI_HEADERIMAGE "${CMAKE_CURRENT_SOURCE_DIR}/docs/img\\\\installer_header.bmp")
  set(CPACK_NSIS_MUI_ICON "${CMAKE_CURRENT_SOURCE_DIR}/docs/img\\\\logo.ico")
  set(CPACK_NSIS_MUI_UNICON "${CMAKE_CURRENT_SOURCE_DIR}/docs/img\\\\logo.ico")
  set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
  set(CPACK_NSIS_MENU_LINKS "bin/multifebe" "MultiFEBE"
						                "doc/manual.pdf" "Manual")
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_GENERATOR NSIS64)
endif (WIN32)

include(CPack)
