# Pawsey Supercomputing Research Centre
#
# Test building executable with profile_util library 

set(tests
    test_profile_util
)

set(ctests
    test_profile_util_c_api
)

foreach(test ${tests})
  add_executable(${test} ${test}.cpp)
  if (PU_ENABLE_HIP)
    set_source_files_properties(${test}.cpp PROPERTIES LANGUAGE HIP)
  endif()

  target_link_libraries(${test} profile_util ${PU_LIBS})
  if (PU_LINK_FLAGS)
    set_target_properties(${test} PROPERTIES LINK_FLAGS ${PU_LINK_FLAGS})
  endif()
endforeach()

if (PU_ENABLE_C_API)
  foreach(test ${ctests})
    message(STATUS "Building C test: ${test}")
    set_source_files_properties(${test}.c PROPERTIES LANGUAGE C)
    add_executable(${test} ${test}.c)
    # if (PU_ENABLE_HIP)
    #   set_source_files_properties(${test}.c PROPERTIES LANGUAGE HIP)
    # endif()

    target_link_libraries(${test} profile_util ${PU_LIBS})
    if (PU_LINK_FLAGS)
      set_target_properties(${test} PROPERTIES LINK_FLAGS ${PU_LINK_FLAGS})
    endif()
  endforeach()
endif()