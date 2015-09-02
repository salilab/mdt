# - Try to find Modeller
# Once done, this will define
#
#  MODELLER_FOUND - system has Modeller
#  MODELLER_INCLUDE_DIRS - the Modeller include directories
#  MODELLER_LIBRARIES - link these to use Modeller

find_package(PkgConfig)
pkg_check_modules(PC_MODELLER QUIET modeller)

find_library(MODELLER_LIBRARIES
    NAMES modeller
    HINTS ${PC_MODELLER_LIBDIR}
          ${PC_MODELLER_LIBRARY_DIRS}
)

find_path(MODELLER_INCLUDE_DIR
    NAMES modeller.h
    HINTS ${PC_MODELLER_INCLUDEDIR}
          ${PC_MODELLER_INCLUDE_DIRS}
    PATH_SUFFIXES modeller
)

find_path(MODELLER_ARCH_INCLUDE_DIR
    NAMES fortran-pointer-types.h
    HINTS ${PC_MODELLER_INCLUDEDIR}
          ${PC_MODELLER_INCLUDE_DIRS}
    PATH_SUFFIXES modeller
)
set(MODELLER_INCLUDE_DIRS ${MODELLER_INCLUDE_DIR} ${MODELLER_ARCH_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MODELLER REQUIRED_VARS MODELLER_INCLUDE_DIRS
                MODELLER_LIBRARIES ${ADDITIONAL_REQUIRED_VARS})
