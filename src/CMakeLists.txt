# set the source dir

set(PU_SOURCES
    "${git_revision_cpp}"
    mem_util.cpp
    thread_affinity_util.cpp
    time_util.cpp
    profile_util.cpp
)

if (PU_ENABLE_SHARED_LIB)
    add_library(profile_util SHARED ${PU_SOURCES})
else()
    add_library(profile_util STATIC ${PU_SOURCES})
endif()
target_compile_definitions(profile_util PUBLIC ${PU_DEFINES})
if (PU_CXX_FLAGS)
    target_compile_options(profile_util PUBLIC ${PU_CXX_FLAGS})
endif()
set_target_properties(profile_util PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (BUILD_TESTING AND PU_ENABLE_TESTS)
	add_subdirectory(tests)
endif()

install(TARGETS profile_util 
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if (PU_ENABLE_PYTHON_INTERFACE)
    pybind11_add_module(py_profile_util
    profile_util_pyinterface.cpp 
    profile_util.cpp 
    mem_util.cpp 
    thread_affinity_util.cpp
    time_util.cpp
    "${git_revision_cpp}")
    #install(TARGETS profile_util DESTINATION .)
endif()
if (PU_ENABLE_HIP)
    set_source_files_properties(mem_util.cpp PROPERTIES LANGUAGE HIP)
    set_source_files_properties(thread_affinity_util.cpp PROPERTIES LANGUAGE HIP)
    set_source_files_properties(time_util.cpp PROPERTIES LANGUAGE HIP)
    set_source_files_properties(profile_util.cpp PROPERTIES LANGUAGE HIP)
endif()
