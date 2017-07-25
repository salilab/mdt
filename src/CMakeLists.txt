include(GenerateExportHeader)
file(GLOB_RECURSE cfiles "*.c")

configure_file("mdt_version.c.in"
               "${PROJECT_BINARY_DIR}/src/mdt_version.c")
list(APPEND cfiles "${PROJECT_BINARY_DIR}/src/mdt_version.c")

include_directories(BEFORE "${CMAKE_BINARY_DIR}/src")
include_directories(".")
add_library(mdt-lib SHARED ${cfiles})
set_target_properties(mdt-lib PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      OUTPUT_NAME "mdt"
                      C_VISIBILITY_PRESET hidden)
set_property(TARGET mdt-lib PROPERTY COMPILE_DEFINITIONS
             "MDT_SIZEOF_POINTER=${CMAKE_SIZEOF_VOID_P}")
generate_export_header(mdt-lib
                       EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/src/mdt_config.h"
                       EXPORT_MACRO_NAME MDTDLLEXPORT
                       NO_EXPORT_MACRO_NAME MDTDLLLOCAL)

target_link_libraries(mdt-lib
                      ${EXTRA_LIBS})
INSTALL(TARGETS mdt-lib DESTINATION ${CMAKE_INSTALL_LIBDIR})
