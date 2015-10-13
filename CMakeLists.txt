cmake_minimum_required(VERSION 2.8)
project(MDT)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)

set(MDT_VERSION "git")

# Use @rpath on Mac by default
if(POLICY CMP0042)
cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

enable_testing()

include(GNUInstallDirs)
if (NOT DEFINED CMAKE_INSTALL_PYTHONDIR)
  set(CMAKE_INSTALL_PYTHONDIR "${CMAKE_INSTALL_LIBDIR}/python/site_packages")
endif()

find_package(HDF5 1.8 REQUIRED COMPONENTS C HL)
set(EXTRA_LIBS ${EXTRA_LIBS} ${HDF5_LIBRARIES})
include_directories(${HDF5_INCLUDE_DIRS})

find_package(GLIB REQUIRED)
set(EXTRA_LIBS ${EXTRA_LIBS} ${GLIB_LIBRARIES})
include_directories(${GLIB_INCLUDE_DIRS})

find_package(MODELLER REQUIRED)
set(EXTRA_LIBS ${EXTRA_LIBS} ${MODELLER_LIBRARIES})
include_directories(${MODELLER_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(pyext)
add_subdirectory(data)
add_subdirectory(test)
