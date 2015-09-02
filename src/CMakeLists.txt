file(GLOB_RECURSE cfiles "${PROJECT_SOURCE_DIR}/src/*.c")

configure_file("${PROJECT_SOURCE_DIR}/src/mdt_version.c.in"
               "${PROJECT_BINARY_DIR}/src/mdt_version.c")
list(APPEND cfiles "${PROJECT_BINARY_DIR}/src/mdt_version.c")

include_directories("${PROJECT_SOURCE_DIR}/src/")
add_library(mdt-lib SHARED ${cfiles})
set_target_properties(mdt-lib PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      OUTPUT_NAME "mdt")
set_property(TARGET mdt-lib PROPERTY COMPILE_DEFINITIONS
             "MDT_EXPORTS;MDT_SIZEOF_POINTER=${CMAKE_SIZEOF_VOID_P}")
target_link_libraries(mdt-lib
                      ${EXTRA_LIBS})
INSTALL(TARGETS mdt-lib DESTINATION ${CMAKE_INSTALL_LIBDIR})
