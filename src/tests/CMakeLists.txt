# Pawsey Supercomputing Research Centre
#
# Test building executable with profile_util library 

set(tests
    test_profile_util
    test_affinity
)
set(gputests
    test_gpu
)
set(mpitests
    test_mpi_comm
    test_mpi_io
    test_mpi_compute
)
set(gpumpitests
    test_gpu_mpi_comm
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

if (PU_ENABLE_MPI)
  foreach(test ${mpitests})
    add_executable(${test} ${test}.cpp)
    if (PU_ENABLE_HIP)
      set_source_files_properties(${test}.cpp PROPERTIES LANGUAGE HIP)
    endif()
    target_link_libraries(${test} profile_util ${PU_LIBS})
    if (PU_LINK_FLAGS)
      set_target_properties(${test} PROPERTIES LINK_FLAGS ${PU_LINK_FLAGS})
    endif()
  endforeach()
endif()

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

if (PU_ENABLE_HIP OR PU_ENABLE_CUDA)
  foreach(test ${gputests})
    add_executable(${test} ${test}.cpp)
    if (PU_ENABLE_HIP)
      set_source_files_properties(${test}.cpp PROPERTIES LANGUAGE HIP)
    endif()
    target_link_libraries(${test} profile_util ${PU_LIBS})
    if (PU_LINK_FLAGS)
      set_target_properties(${test} PROPERTIES LINK_FLAGS ${PU_LINK_FLAGS})
    endif()
  endforeach()
  target_sources(test_gpu PRIVATE test_gpu_kernels.cpp)
  if (PU_ENABLE_HIP)
    set_source_files_properties(test_gpu_kernels.cpp PROPERTIES LANGUAGE HIP)
  elseif (NOT CMAKE_CXX_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
    set_source_files_properties(test_gpu_kernels.cpp PROPERTIES LANGUAGE CUDA)
    set_target_properties(test_gpu PROPERTIES LINKER_LANGUAGE CXX CUDA_RUNTIME_LIBRARY Shared)
  endif()
endif()

if (PU_ENABLE_MPI AND (PU_ENABLE_HIP OR PU_ENABLE_CUDA))
  foreach(test ${gpumpitests})
    add_executable(${test} ${test}.cpp)
    if (PU_ENABLE_HIP)
      set_source_files_properties(${test}.cpp PROPERTIES LANGUAGE HIP)
    endif()
    target_link_libraries(${test} profile_util ${PU_LIBS})
    if (PU_LINK_FLAGS)
      set_target_properties(${test} PROPERTIES LINK_FLAGS ${PU_LINK_FLAGS})
    endif()
  endforeach()
endif()
