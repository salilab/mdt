find_package(SWIG 1.3.40 REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonInterp 2.4 REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
include_directories(BEFORE "${PROJECT_SOURCE_DIR}/src")

FILE(GLOB SWIG_INCLUDES "${PROJECT_SOURCE_DIR}/swig/*.i")

SET(CMAKE_SWIG_FLAGS "-nodefaultctor" "-nodefaultdtor")
SET(SWIG_MODULE_mdt_EXTRA_DEPS ${SWIG_INCLUDES})

SWIG_ADD_MODULE(mdt PYTHON mdt.i)

if(APPLE)
  set_target_properties("${SWIG_MODULE_mdt_REAL_NAME}" PROPERTIES LINK_FLAGS "-flat_namespace -undefined suppress")
endif(APPLE)

# Don't link in libpython (that prevents the extension from being used in
# later compatible versions of Python), except on Windows where all DLL
# symbols must be defined
IF(WIN32)
  SWIG_LINK_LIBRARIES(mdt ${PYTHON_LIBRARIES} mdt-lib)
ELSE(WIN32)
  SWIG_LINK_LIBRARIES(mdt mdt-lib)
ENDIF(WIN32)

set_target_properties("${SWIG_MODULE_mdt_REAL_NAME}"
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

INSTALL(TARGETS ${SWIG_MODULE_mdt_REAL_NAME} DESTINATION ${CMAKE_INSTALL_PYTHONDIR})
INSTALL(FILES "${CMAKE_SOURCE_DIR}/pyext/mdt/__init__.py" DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/mdt/")
INSTALL(FILES "${CMAKE_SOURCE_DIR}/pyext/mdt/features.py" DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/mdt/")

set(mdt_SWIG _mdt CACHE INTERNAL "" FORCE)
