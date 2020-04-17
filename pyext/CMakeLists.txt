find_package(SWIG 1.3.40 REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonInterp 2.4 REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)

find_package(NumPy)
if(NUMPY_FOUND)
  include_directories(${NUMPY_INCLUDE_DIRS})
  add_definitions("-DMDT_WITH_NUMPY")
endif()

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
include_directories(BEFORE "${PROJECT_SOURCE_DIR}/src")
include_directories(BEFORE "${CMAKE_BINARY_DIR}/src")

FILE(GLOB SWIG_INCLUDES "*.i")

SET(CMAKE_SWIG_FLAGS "-nodefaultctor" "-nodefaultdtor" "-noproxy")
SET(SWIG_MODULE__mdt_EXTRA_DEPS ${SWIG_INCLUDES})

# cmake 3.8 replaced SWIG_ADD_MODULE with SWIG_ADD_LIBRARY
IF (COMMAND SWIG_ADD_LIBRARY)
  SWIG_ADD_LIBRARY(_mdt LANGUAGE python TYPE MODULE SOURCES mdt.i)
ELSE()
  SWIG_ADD_MODULE(_mdt PYTHON mdt.i)
ENDIF()

if(APPLE)
  set_target_properties("${SWIG_MODULE__mdt_REAL_NAME}" PROPERTIES LINK_FLAGS "-flat_namespace -undefined suppress")
endif(APPLE)

# Don't link in libpython (that prevents the extension from being used in
# later compatible versions of Python), except on Windows where all DLL
# symbols must be defined
IF(WIN32)
  SWIG_LINK_LIBRARIES(_mdt ${PYTHON_LIBRARIES} mdt-lib)
ELSE(WIN32)
  SWIG_LINK_LIBRARIES(_mdt mdt-lib)
ENDIF(WIN32)

set_target_properties("${SWIG_MODULE__mdt_REAL_NAME}"
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

INSTALL(TARGETS ${SWIG_MODULE__mdt_REAL_NAME} DESTINATION ${CMAKE_INSTALL_PYTHONDIR})
INSTALL(FILES "mdt/__init__.py" DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/mdt/")
INSTALL(FILES "mdt/features.py" DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/mdt/")
