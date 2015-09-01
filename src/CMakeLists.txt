file(GLOB_RECURSE cfiles "${PROJECT_SOURCE_DIR}/src/*.c")

add_library(mdt-lib SHARED ${cfiles})
set_target_properties(mdt-lib PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      OUTPUT_NAME "mdt")
set_property(TARGET mdt-lib PROPERTY COMPILE_DEFINITIONS
             "MDT_EXPORTS;MDT_SIZEOF_POINTER=${CMAKE_SIZEOF_VOID_P}")
target_link_libraries(mdt-lib
                      ${EXTRA_LIBS})
INSTALL(TARGETS mdt-lib DESTINATION ${CMAKE_INSTALL_LIBDIR})
