#
# How we find MPI and set it up
#
macro(pu_find_mpi)
    message("MPI enabled, finding MPI ... ")
    find_package(MPI)
    if (MPI_FOUND)
        list(APPEND PU_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
        list(APPEND PU_LIBS ${MPI_CXX_LIBRARIES})
        list(APPEND PU_CXX_FLAGS ${MPI_CXX_FLAGS})
        list(APPEND PU_LINK_FLAGS ${MPI_CXX_FLAGS})
        list(APPEND PU_DEFINES "_MPI")
        set(PU_HAS_MPI Yes)
    else()
        message(SEND_ERROR "MPI enabled but not found. Please check configuration or disable MPI")
    endif()
endmacro()

macro(pu_find_openmp)
    message("OpenMP enabled, finding OpenMP ... ")
    find_package(OpenMP)
    #FindOpenMP()
    if (OpenMP_FOUND)
        list(APPEND PU_CXX_FLAGS ${OpenMP_CXX_FLAGS})
        list(APPEND PU_LINK_FLAGS ${OpenMP_CXX_FLAGS})
        #list(APPEND PU_DEFINES ${_OPENMP})
        set(PU_HAS_OPENMP Yes)
    else()
        message(SEND_ERROR "OpenMP enabled but not found. Please check configuration or disable OpenMP")
    endif()
endmacro()

macro(pu_find_hip)
    # Find hip
    message("HIP enabled, finding HIP ... ")
    find_package(HIP)
    if (HIP_FOUND)
        # Link with HIP
        enable_language(HIP)
        if(NOT DEFINED CMAKE_HIP_STANDARD)
            set(CMAKE_HIP_STANDARD 17)
            set(CMAKE_HIP_STANDARD_REQUIRED ON)
        endif()
        list(APPEND PU_DEFINES "_HIP")
        set(PU_HAS_HIP Yes)
        add_compile_options("-fPIE")
        target_link_libraries(hip::device)
        if (ENABLE_HIP_AMD)
            message("Using AMD HIP ")
            add_definitions("-D__HIP_PLATFORM_AMD__")
            set(GPU_TARGETS "gfx90a" CACHE STRING "GPU targets to compile for")
        endif()
    else()
        message(SEND_ERROR "HIP enabled but not found. Please check configuration or disable HIP")
    endif()
endmacro()

macro(pu_find_cuda)
    # Find hip
    message(STATUS "CUDA enabled, finding CUDA ... ")
    find_package(CUDA)
    if (CUDA_FOUND)
        enable_language(CUDA)
        list(APPEND PU_DEFINES "_CUDA")
        set(PU_HAS_CUDA Yes)
        if(NOT DEFINED CMAKE_CUDA_STANDARD)
            set(CMAKE_CUDA_STANDARD 17)
            set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        endif()
    else()
        message(SEND_ERROR "CUDA enabled but not found. Please check configuration or disable CUDA")
    endif()
endmacro()

macro(pu_mpi)
    set(PU_HAS_MPI No)
    if (PU_ENABLE_MPI)
        pu_find_mpi()
    endif()
endmacro()

macro(pu_openmp)
    set(PU_HAS_OPENMP No)
    if (PU_ENABLE_OPENMP)
        pu_find_openmp()
    endif()
endmacro()

macro(pu_hip)
    set(PU_HAS_HIP No)
    if (PU_ENABLE_HIP)
        pu_find_hip()
    endif()
endmacro()

macro(pu_cuda)
    set(PU_HAS_CUDA No)
    if (PU_ENABLE_CUDA)
        pu_find_cuda()
    endif()
endmacro()

macro(pu_report feature)

    # Output feature name and underscore it in the next line
    message("\n${feature}")
    string(REGEX REPLACE "." "-" _underscores ${feature})
    message("${_underscores}\n")

    set(_args "${ARGN}")
    list(LENGTH _args _nargs)
    math(EXPR _nargs "${_nargs} - 1")
    foreach(_idx RANGE 0 ${_nargs} 2)
        # Items in the list come with a message first, then the variable name
        list(GET _args ${_idx} _msg)
        math(EXPR _idx2 "${_idx} + 1")
        list(GET _args ${_idx2} _varname)

        # We try to keep things up to 80 cols
        string(LENGTH ${_msg} _len)
        math(EXPR _nspaces "75 - ${_len}")
        string(RANDOM LENGTH ${_nspaces} _spaces)
        string(REGEX REPLACE "." " " _spaces "${_spaces}")
        string(CONCAT _msg "${_msg}" ${_spaces})
        message(" ${_msg} ${VR_HAS_${_varname}}")
    endforeach()
endmacro()

