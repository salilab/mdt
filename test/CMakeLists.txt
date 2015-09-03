file(GLOB_RECURSE tests "test_*.py")

foreach(test ${tests})
  get_filename_component(name ${test} NAME)
  add_test(NAME ${name} COMMAND python ${test}
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  set_tests_properties(${name} PROPERTIES TIMEOUT 5)
  if (DEFINED ENV{PYTHONPATH})
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pyext:${CMAKE_SOURCE_DIR}/pyext:$ENV{PYTHONPATH}")
  else()
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pyext:${CMAKE_SOURCE_DIR}/pyext")
  endif()
endforeach()
